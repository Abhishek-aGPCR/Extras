/* SYNE1 isoform browser - variant engine.
 *
 * One implementation shared by the HTML browser (window.VariantEngine) and the command-line
 * report (scripts/07_variant_report.js, Node). It works on the browser's data model: a "view"
 * is one genomic locus ({key, species, symbol, assembly, chrom, chrom_acc, region, seq (+ strand),
 * giant_ref, domains, transcripts}), a transcript has exons in transcript order, strand and a
 * genomic CDS range that includes the stop codon.
 *
 * What it does
 *   parse()      c./n./g./chr:pos/VCF/p./exon-relative descriptions
 *   resolve()    -> genomic variant(s) on the + strand, with a reference-allele check
 *   analyse()    -> per-isoform location, HGVS, codon/amino-acid change, splice effects,
 *                   premature termination codon (PTC) and NMD prediction, protein domains lost,
 *                   orthologous positions in the other species
 *   figureSVG()  -> gene / mRNA-fate / protein cartoon as a standalone SVG
 *
 * NMD rules (Lindeboom, Supek & Lehner 2016, Nat Genet 48:1112; Nagy & Maquat 1998): a PTC more
 * than 50-55 nt upstream of the last exon-exon junction triggers NMD; PTCs in the last exon or
 * within ~50 nt of the last junction escape; start-proximal PTCs (<150 nt from the start codon)
 * and PTCs in long exons (>407 nt) can escape. Predictions, not measurements.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.VariantEngine = factory();
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const COMP = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };
  const rc = s => { let o = ''; for (let i = s.length - 1; i >= 0; i--) o += COMP[s[i]] || 'N'; return o; };
  const CODON = {};
  (function () { const b = 'TCAG', aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; let k = 0;
    for (const x of b) for (const y of b) for (const z of b) CODON[x + y + z] = aa[k++]; })();
  const translate = nt => { let p = ''; for (let i = 0; i + 2 < nt.length; i += 3) p += CODON[nt.substr(i, 3)] || 'X'; return p; };
  const AA3 = { A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys',
    M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Ter', X: 'Xaa' };
  const AA1 = {}; for (const k in AA3) AA1[AA3[k].toUpperCase()] = k; AA1.STOP = '*';
  const fmt = n => Number(n).toLocaleString('en-US');
  const esc = s => String(s == null ? '' : s).replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));

  /* ------------------------------------------------------------------ geometry */
  function geom(t) {
    if (t._vg) return t._vg;
    const cum = [0]; for (const [s, e] of t.exons) cum.push(cum[cum.length - 1] + e - s + 1);
    const g = { cum, len: cum[cum.length - 1], cds0: -1, cds1: -1 };
    t._vg = g;
    if (t.cds) {
      const a = t.strand === 1 ? t.cds[0] : t.cds[1], b = t.strand === 1 ? t.cds[1] : t.cds[0];
      g.cds0 = g2i(t, a); g.cds1 = g2i(t, b);
    }
    const first = t.exons[0], last = t.exons[t.exons.length - 1];
    g.tss = t.strand === 1 ? first[0] : first[1];
    g.tes = t.strand === 1 ? last[1] : last[0];
    g.lo = Math.min(...t.exons.map(x => x[0])); g.hi = Math.max(...t.exons.map(x => x[1]));
    return g;
  }
  function exonOf(t, g) { for (let k = 0; k < t.exons.length; k++) { const [s, e] = t.exons[k]; if (g >= s && g <= e) return k; } return -1; }
  function g2i(t, g) { const cum = t._vg ? t._vg.cum : null; let acc = 0;
    for (let k = 0; k < t.exons.length; k++) { const [s, e] = t.exons[k]; if (g >= s && g <= e) return (cum ? cum[k] : acc) + (t.strand === 1 ? g - s : e - g); acc += e - s + 1; }
    return -1; }
  function i2g(t, i) { const cum = geom(t).cum; if (i < 0 || i >= cum[cum.length - 1]) return null;
    let k = 0; while (cum[k + 1] <= i) k++; const [s, e] = t.exons[k]; const off = i - cum[k]; return t.strand === 1 ? s + off : e - off; }
  function exonOfIdx(t, i) { const cum = geom(t).cum; let k = 0; while (k < t.exons.length - 1 && cum[k + 1] <= i) k++; return k; }
  function gsub(v, a, b, strand) { const s = v.seq.slice(a - v.region.start, b - v.region.start + 1); return strand === 1 ? s : rc(s); }
  const MRNA = new WeakMap();
  function mrna(v, t) { if (!MRNA.has(t)) MRNA.set(t, t.exons.map(([s, e]) => gsub(v, s, e, t.strand)).join('')); return MRNA.get(t); }
  function byId(v, id) { if (!id) return null; for (const t of v.transcripts) if (t.id === id || t.id.split('.')[0] === String(id).split('.')[0]) return t; return null; }

  function cpos(t, i) {
    const G = geom(t); if (G.cds0 < 0) return String(i + 1);
    if (i < G.cds0) return '-' + (G.cds0 - i);
    if (i > G.cds1) return '*' + (i - G.cds1);
    return String(i - G.cds0 + 1);
  }
  /** where genomic position g falls in transcript t */
  function locate(t, g) {
    const G = geom(t), pre = G.cds0 >= 0 ? 'c.' : 'n.';
    const k = exonOf(t, g);
    if (k >= 0) { const i = g2i(t, g); return { type: 'exon', exon: k, i, pos: pre + cpos(t, i) }; }
    const up = t.strand === 1 ? g < G.tss : g > G.tss, down = t.strand === 1 ? g > G.tes : g < G.tes;
    if (up) { const d = Math.abs(g - G.tss); return { type: 'upstream', d, pos: pre + (G.cds0 >= 0 ? '-' + (G.cds0 + d) : '-' + d) }; }
    if (down) { const d = Math.abs(g - G.tes); return { type: 'downstream', d, pos: pre + (G.cds0 >= 0 ? '*' + (G.len - 1 - G.cds1 + d) : '*' + d) }; }
    for (let k2 = 0; k2 < t.exons.length - 1; k2++) {
      const A = t.exons[k2], B = t.exons[k2 + 1];
      const lo = t.strand === 1 ? A[1] : B[1], hi = t.strand === 1 ? B[0] : A[0];
      if (g > lo && g < hi) {
        const upEnd = t.strand === 1 ? A[1] : A[0], dnStart = t.strand === 1 ? B[0] : B[1];
        const dUp = Math.abs(g - upEnd), dDown = Math.abs(dnStart - g);
        const pos = dUp <= dDown ? pre + cpos(t, g2i(t, upEnd)) + '+' + dUp : pre + cpos(t, g2i(t, dnStart)) + '-' + dDown;
        return { type: 'intron', intron: k2, dUp, dDown, len: hi - lo - 1, pos };
      }
    }
    return { type: 'unknown', pos: '' };
  }

  /* ------------------------------------------------------------------ parsing */
  const NT = '[ACGTNacgtn]';
  function parse(text) {
    let s = String(text || '').trim().replace(/ /g, ' ');
    if (!s) return { error: 'Enter a variant, for example c.23560G>T, chr6:152176461C>A or p.Glu7854*.' };
    let acc = null, m;
    if ((m = s.match(/^([A-Za-z]{2}_\d+(?:\.\d+)?|ENS[A-Z]*[TG]\d+(?:\.\d+)?|LRG_\d+(?:t\d+)?)\s*(?:\([^)]*\))?\s*:\s*/))) { acc = m[1]; s = s.slice(m[0].length); }
    s = s.replace(/^\((.*)\)$/, '$1').trim();
    // c. / n.
    if ((m = s.match(new RegExp('^([cn])\\.\\s*\\(?([-*]?\\d+)([+-]\\d+)?(?:_([-*]?\\d+)([+-]\\d+)?)?\\)?\\s*(.*)$', 'i')))) {
      const edit = parseEdit(m[6]); if (edit.error) return edit;
      return Object.assign({ kind: m[1].toLowerCase(), acc, p1: m[2], o1: m[3] ? +m[3] : 0, p2: m[4] || null, o2: m[5] ? +m[5] : 0 }, edit);
    }
    // g.
    if ((m = s.match(/^g\.\s*(\d+)(?:_(\d+))?\s*(.*)$/i))) {
      const edit = parseEdit(m[3]); if (edit.error) return edit;
      return Object.assign({ kind: 'g', acc, start: +m[1], end: m[2] ? +m[2] : +m[1] }, edit);
    }
    // p.
    if ((m = s.match(/^p\.\s*\(?\s*([A-Za-z]{3}|[A-Za-z*])\s*(\d+)\s*([A-Za-z]{3}|[A-Za-z*=]|fs\S*)\s*\)?$/i)) ||
        (m = s.match(/^([A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWY])(\d+)([A-Z][a-z]{2}|[ACDEFGHIKLMNPQRSTVWYX*]|=)$/))) {
      const a1 = aaCode(m[1]), a2 = /^fs/i.test(m[3]) ? 'fs' : aaCode(m[3]);
      if (!a1 || !a2) return { error: `Unrecognised amino acid in "${s}". Use one- or three-letter codes (Glu or E; Ter, * or X for a stop).` };
      return { kind: 'p', acc, ref: a1, pos: +m[2], alt: a2 };
    }
    // exon-relative: exon 130:99G>T
    if ((m = s.match(new RegExp('^(?:exon|ex|e)\\s*(\\d+)\\s*[:@ ,]\\s*(?:pos(?:ition)?\\s*)?(\\d+)\\s*(?:(' + NT + ')\\s*>\\s*(' + NT + '))?$', 'i')))) {
      return { kind: 'exon', acc, exon: +m[1], pos: +m[2], ref: m[3] ? m[3].toUpperCase() : null, alt: m[4] ? m[4].toUpperCase() : null, type: m[4] ? 'sub' : 'pos' };
    }
    // genomic chr:pos ref>alt, VCF-like 6-152176461-C-A, "chr6 152176461 C A"
    if ((m = s.match(new RegExp('^(?:chr)?([0-9]{1,2}|[XYM]|MT)\\s*[:\\-\\s]\\s*([\\d,]+)\\s*(?:[:\\-\\s]?\\s*(' + NT + '+)\\s*(?:>|[:\\-\\s/])\\s*(' + NT + '+|-))?$', 'i')))) {
      const pos = +m[2].replace(/,/g, '');
      if (!m[3]) return { kind: 'g', acc, chrom: m[1], start: pos, end: pos, type: 'pos' };
      const ref = m[3].toUpperCase(), alt = m[4] === '-' ? '' : m[4].toUpperCase();
      return { kind: 'g', acc, chrom: m[1], start: pos, end: pos + ref.length - 1, type: ref.length === 1 && alt.length === 1 ? 'sub' : 'vcf', ref, alt };
    }
    return { error: `Could not read "${text}". Examples: c.23560G>T, NM_182961.4:c.23560G>T, c.123+1G>A, c.100del, n.24125G>T, chr6:152176461C>A, 6-152176461-C-A, g.152176461C>A, p.Glu7854*, E7854X, exon 130:99G>T.` };
  }
  function parseEdit(e) {
    e = (e || '').trim(); let m;
    if (!e) return { type: 'pos' };
    if ((m = e.match(new RegExp('^(' + NT + ')\\s*>\\s*(' + NT + ')$')))) return { type: 'sub', ref: m[1].toUpperCase(), alt: m[2].toUpperCase() };
    if ((m = e.match(new RegExp('^delins(' + NT + '+)$', 'i')))) return { type: 'delins', alt: m[1].toUpperCase() };
    if ((m = e.match(new RegExp('^del(' + NT + '*)$', 'i')))) return { type: 'del', delseq: m[1].toUpperCase() };
    if ((m = e.match(new RegExp('^dup(' + NT + '*)$', 'i')))) return { type: 'dup', dupseq: m[1].toUpperCase() };
    if ((m = e.match(new RegExp('^ins(' + NT + '+)$', 'i')))) return { type: 'ins', ins: m[1].toUpperCase() };
    if (e === '=') return { type: 'pos' };
    return { error: `Could not read the change "${e}". Supported: N>N, del, dup, insN, delinsN.` };
  }
  function aaCode(x) { if (!x) return null; if (x === '=' ) return '='; if (x === '*' || /^x$/i.test(x)) return '*';
    if (x.length === 1) return /[ACDEFGHIKLMNPQRSTVWY]/i.test(x) ? x.toUpperCase() : null;
    return AA1[x.toUpperCase()] || null; }

  /* ------------------------------------------------------------------ resolve to + strand genomic */
  function cToIdx(t, p) {
    const G = geom(t);
    if (G.cds0 < 0) return null;
    if (p[0] === '*') return G.cds1 + (+p.slice(1));
    const n = +p; if (n < 0) return G.cds0 + n; return G.cds0 + n - 1;
  }
  function idxToG(t, i, off) {
    // positions beyond the transcript ends are extrapolated along the genome
    const G = geom(t); let g;
    if (i < 0) g = G.tss + i * t.strand;
    else if (i >= G.len) g = G.tes + (i - G.len + 1) * t.strand;
    else g = i2g(t, i);
    return g + (off || 0) * t.strand;
  }
  function resolve(v, P, txId) {
    if (P.error) return { error: P.error };
    const notes = [];
    const pick = () => { const t = byId(v, P.acc) || byId(v, txId) || byId(v, v.giant_ref); if (P.acc && !byId(v, P.acc)) notes.push(`${P.acc} is not among the ${v.symbol} isoforms of this view; ${t.id} was used instead.`); return t; };
    const plus = (b, strand) => strand === 1 ? b : rc(b);
    const base = g => (g >= v.region.start && g <= v.region.end) ? v.seq[g - v.region.start] : null;
    const range = (a, b) => (a >= v.region.start && b <= v.region.end) ? v.seq.slice(a - v.region.start, b - v.region.start + 1) : null;
    let out = [];
    if (P.kind === 'c' || P.kind === 'n') {
      const t = pick();
      const G = geom(t);
      const idxOf = p => P.kind === 'n' ? +p - 1 : cToIdx(t, p);
      for (const [p, o] of [[P.p1, P.o1], [P.p2, P.o2]]) {
        if (p == null || !o) continue;
        const i = idxOf(p); if (i == null) break;
        const k = exonOfIdx(t, i), first = G.cum[k], last = G.cum[k + 1] - 1;
        if ((o > 0 && i !== last) || (o < 0 && i !== first)) {
          const fC = cpos(t, first), lC = cpos(t, last);
          return { error: `${P.kind}.${p}${o > 0 ? '+' : ''}${o} is not valid HGVS: an intronic offset must start from an exon boundary. ${P.kind}.${p} lies inside exon ${k + 1} of ${t.id}, which spans ${P.kind}.${fC} to ${P.kind}.${lC}; use ${P.kind}.${lC}+N for the intron after it or ${P.kind}.${fC}-N for the intron before it.` };
        }
      }
      const toG = (p, o) => { const i = idxOf(p); if (i == null) return null; return idxToG(t, i, o); };
      const g1 = toG(P.p1, P.o1); if (g1 == null) return { error: `${t.id} has no CDS; use n. (transcript) coordinates.` };
      const g2 = P.p2 ? toG(P.p2, P.o2) : g1;
      const a = Math.min(g1, g2), b = Math.max(g1, g2);
      if (range(a, b) == null) return { error: 'The position lies outside the downloaded locus.' };
      const S = t.strand;
      if (P.type === 'sub' || P.type === 'pos') {
        const refP = base(a);
        const v1 = { start: a, end: a, ref: refP, alt: P.type === 'sub' ? plus(P.alt, S) : null, type: P.type, tx: t.id,
          given_ref: P.ref ? plus(P.ref, S) : null };
        out.push(v1);
      } else if (P.type === 'del') {
        out.push({ start: a, end: b, ref: range(a, b), alt: '', type: 'del', tx: t.id, given_ref: P.delseq ? plus(P.delseq, S) : null });
      } else if (P.type === 'delins') {
        out.push({ start: a, end: b, ref: range(a, b), alt: plus(P.alt, S), type: 'delins', tx: t.id });
      } else if (P.type === 'dup') {
        const seg = range(a, b); out.push({ start: b, end: b, ref: base(b), alt: base(b) + seg, type: 'dup', tx: t.id, given_ref: P.dupseq ? plus(P.dupseq, S) : null, display_ref: seg, dupStart: a, dupEnd: b });
      } else if (P.type === 'ins') {
        if (b - a !== 1) return { error: 'An insertion needs two adjacent flanking positions, e.g. c.123_124insA.' };
        out.push({ start: a, end: a, ref: base(a), alt: base(a) + plus(P.ins, S), type: 'ins', tx: t.id });
      }
    } else if (P.kind === 'g') {
      if (P.chrom && String(P.chrom).replace(/^chr/i, '') !== String(v.chrom)) notes.push(`Chromosome ${P.chrom} given; this locus is on chr${v.chrom}.`);
      if (P.acc && v.chrom_acc && P.acc.split('.')[0] !== v.chrom_acc.split('.')[0] && /^NC_/.test(P.acc)) notes.push(`${P.acc} differs from this view's ${v.chrom_acc}.`);
      const a = P.start, b = P.end;
      if (range(a, b) == null) return { error: `Position ${fmt(a)} lies outside the downloaded ${v.symbol} locus (chr${v.chrom}:${fmt(v.region.start)}-${fmt(v.region.end)}, ${v.assembly}).` };
      if (P.type === 'sub' || P.type === 'pos') out.push({ start: a, end: a, ref: base(a), alt: P.alt || null, type: P.type, given_ref: P.ref || null });
      else if (P.type === 'vcf') {
        // left-anchored VCF style: trim shared prefix
        let ref = P.ref, alt = P.alt, s = a;
        while (ref.length && alt.length && ref[0] === alt[0] && (ref.length > 1 || alt.length > 1)) { ref = ref.slice(1); alt = alt.slice(1); s++; }
        if (!ref.length) out.push({ start: s - 1, end: s - 1, ref: base(s - 1), alt: base(s - 1) + alt, type: 'ins', given_ref: null });
        else out.push({ start: s, end: s + ref.length - 1, ref: range(s, s + ref.length - 1), alt, type: alt.length ? 'delins' : 'del', given_ref: ref });
      } else if (P.type === 'del') out.push({ start: a, end: b, ref: range(a, b), alt: '', type: 'del', given_ref: P.delseq || null });
      else if (P.type === 'delins') out.push({ start: a, end: b, ref: range(a, b), alt: P.alt, type: 'delins' });
      else if (P.type === 'dup') out.push({ start: b, end: b, ref: base(b), alt: base(b) + range(a, b), type: 'dup', dupStart: a, dupEnd: b, display_ref: range(a, b), given_ref: P.dupseq || null });
      else if (P.type === 'ins') out.push({ start: a, end: a, ref: base(a), alt: base(a) + P.ins, type: 'ins' });
    } else if (P.kind === 'p') {
      const t = pick(); const G = geom(t);
      if (G.cds0 < 0) return { error: `${t.id} is non-coding.` };
      const prot = translate(mrna(v, t).slice(G.cds0, G.cds1 + 1));
      if (P.pos < 1 || P.pos > prot.length) return { error: `Residue ${P.pos} is outside ${t.id} (${prot.length - 1} aa).` };
      const refAA = prot[P.pos - 1];
      if (P.ref !== refAA) notes.push(`${t.id} has ${AA3[refAA]}${P.pos}, not ${AA3[P.ref] || P.ref}${P.pos}. Check the isoform: residue numbers differ between isoforms.`);
      if (P.alt === 'fs') return { error: 'A frameshift needs a nucleotide description (e.g. c.23560del): p.(…fs) does not say which bases changed.' };
      const ci = G.cds0 + 3 * (P.pos - 1); const codon = mrna(v, t).substr(ci, 3);
      const cand = [];
      for (let j = 0; j < 3; j++) for (const b of 'ACGT') {
        if (b === codon[j]) continue; const mut = codon.slice(0, j) + b + codon.slice(j + 1);
        const aa = CODON[mut]; if (aa === P.alt || (P.alt === '=' && aa === refAA)) {
          const g = i2g(t, ci + j); cand.push({ start: g, end: g, ref: base(g), alt: t.strand === 1 ? b : COMP[b], type: 'sub', tx: t.id, from_protein: `${codon}>${mut}` });
        }
      }
      if (!cand.length) return { error: `No single-nucleotide change turns ${t.id} codon ${P.pos} (${codon}, ${AA3[refAA]}) into ${AA3[P.alt] || P.alt}.`, notes };
      if (cand.length > 1) notes.push(`${cand.length} single-nucleotide changes encode ${AA3[refAA]}${P.pos}${AA3[P.alt] || P.alt} in ${t.id}; each is analysed below.`);
      out = cand;
    } else if (P.kind === 'exon') {
      const t = pick(); if (P.exon < 1 || P.exon > t.exons.length) return { error: `${t.id} has ${t.exons.length} exons.` };
      const [s, e] = t.exons[P.exon - 1]; if (P.pos < 1 || P.pos > e - s + 1) return { error: `Exon ${P.exon} of ${t.id} is ${e - s + 1} bp long.` };
      const g = t.strand === 1 ? s + P.pos - 1 : e - P.pos + 1;
      out.push({ start: g, end: g, ref: base(g), alt: P.alt ? plus(P.alt, t.strand) : null, type: P.alt ? 'sub' : 'pos', tx: t.id, given_ref: P.ref ? plus(P.ref, t.strand) : null });
    }
    for (const x of out) if (x.type === 'sub' && x.alt && x.alt === x.ref)
      return { error: `No change: the alternative allele (${x.alt}, + strand) is the reference base of ${v.assembly} at ${fmt(x.start)}. Check the alleles and the strand.`, notes };
    for (const x of out) {
      x.ref_ok = x.given_ref ? (x.given_ref === (x.type === 'dup' ? x.display_ref : x.ref)) : null;
      if (x.ref_ok === false) notes.push(`Reference check failed: the input says ${x.given_ref} (+ strand) but the ${v.assembly} genome has ${x.type === 'dup' ? x.display_ref : x.ref} at ${fmt(x.start)}.`);
    }
    return { variants: out, notes };
  }

  /* ------------------------------------------------------------------ NMD */
  function nmd(exonLens, ptcI, cds0) {
    const n = exonLens.length, cum = [0]; for (const L of exonLens) cum.push(cum[cum.length - 1] + L);
    let k = 0; while (k < n - 1 && cum[k + 1] <= ptcI) k++;
    const lastEJ = cum[n - 1], dist = lastEJ - (ptcI + 3);
    const r = { exon: k + 1, nExons: n, distToLastEJ: dist, cdsOffset: ptcI - cds0, exonLen: exonLens[k] };
    if (n === 1) return Object.assign(r, { verdict: 'escape', rule: 'single-exon transcript: no exon junction downstream' });
    if (k === n - 1) return Object.assign(r, { verdict: 'escape', rule: 'PTC in the last exon' });
    if (dist <= 50) return Object.assign(r, { verdict: 'escape', rule: `PTC ${dist} nt upstream of the last exon-exon junction (≤50 nt)` });
    if (dist <= 55) return Object.assign(r, { verdict: 'possible escape', rule: `PTC ${dist} nt upstream of the last exon-exon junction (50-55 nt boundary)` });
    const why = [];
    if (ptcI - cds0 < 150) why.push(`start-proximal PTC (${ptcI - cds0} nt after the start codon; re-initiation can rescue)`);
    if (exonLens[k] > 407) why.push(`PTC in a long exon (${exonLens[k]} nt > 407 nt)`);
    if (why.length) return Object.assign(r, { verdict: 'possible escape', rule: why.join('; ') });
    return Object.assign(r, { verdict: 'NMD', rule: `PTC ${fmt(dist)} nt upstream of the last exon-exon junction (exon ${k + 1} of ${n})` });
  }

  /* ------------------------------------------------------------------ per-transcript consequence */
  function domainsFor(v, t) {
    const mrnaCache = tt => mrna(v, tt);
    const giant = byId(v, v.giant_ref); const doms = v.domains || [];
    if (!giant || !doms.length || !t.cds) return [];
    if (t === giant) return doms.map(d => Object.assign({}, d));
    if (!t._vdom) {
      const Gg = geom(giant), Gt = geom(t);
      const map = new Map();
      for (let r = 0; 3 * r + Gg.cds0 + 2 <= Gg.cds1; r++) map.set(i2g(giant, Gg.cds0 + 3 * r), r + 1);
      const rmap = []; for (let r = 0; 3 * r + Gt.cds0 + 2 <= Gt.cds1; r++) rmap.push(map.get(i2g(t, Gt.cds0 + 3 * r)) || 0);
      const dAt = gr => { if (!gr) return null; for (const d of doms) if (d.type !== 'motif' && gr >= d.start && gr <= d.end) return d; return null; };
      const out = []; let cur = null;
      rmap.forEach((gr, i) => {
        const d = dAt(gr); const key = gr ? (d ? d.label : '') : 'unique';
        if (cur && cur.key === key && (key === 'unique' || gr === cur.lastG + 1)) { cur.end = i + 1; cur.lastG = gr; return; }
        cur = { key, start: i + 1, end: i + 1, lastG: gr, label: key === 'unique' ? 'isoform-specific' : (d ? d.label : ''), type: key === 'unique' ? 'unique' : (d ? d.type : 'none') };
        out.push(cur);
      });
      for (const d of doms) if (d.type === 'motif') { const i = rmap.indexOf(d.start); if (i >= 0) out.push({ start: i + 1, end: i + d.end - d.start + 1, label: d.label, type: 'motif' }); }
      // name the segments that have no counterpart in the reference giant protein
      const prot = translate(mrnaCache(t).slice(Gt.cds0, Gt.cds1 + 1));
      for (const d of out) {
        if (d.type !== 'unique') continue;
        const seg = prot.slice(d.start - 1, d.end);
        if (d.start === 1 && t.category === 'alpha2') { d.label = 'α2-specific N-terminus (exon 1α2)'; }
        else if (d.start === 1) { d.label = 'isoform-specific N-terminus'; d.type = 'extra'; }
        else if (d.end >= rmap.length) { d.label = 'isoform-specific C-terminus'; d.type = 'extra'; }
        else if (/[VIMLA]IPE[STA]PE[AV]/.test(seg) && seg.length <= 30) { d.label = 'DV23 (not in the reference giant)'; d.type = 'dv23'; }
        else { d.label = 'not in the reference giant'; d.type = 'extra'; }
      }
      t._vdom = out.filter(d => d.type !== 'none').map(({ key, lastG, ...d }) => d);
    }
    return t._vdom;
  }
  function translateFrom(seq, i0) { let p = ''; for (let i = i0; i + 2 < seq.length; i += 3) { const a = CODON[seq.substr(i, 3)] || 'X'; p += a; if (a === '*') break; } return p; }

  function consequence(v, t, x) {
    const G = geom(t), S = t.strand, res = { id: t.id, name: t.name, category: t.category, source: t.source, tags: t.tags || [], strand: S };
    const loc = locate(t, x.start); res.loc = loc;
    const locEnd = x.end !== x.start ? locate(t, x.end) : loc;
    const refT = x.ref != null ? (S === 1 ? x.ref : rc(x.ref)) : '', altT = x.alt != null ? (S === 1 ? x.alt : rc(x.alt)) : null;
    // HGVS for this transcript
    if (x.type === 'sub' || x.type === 'pos') res.hgvs = loc.pos + (x.type === 'sub' ? refT + '>' + altT : '');
    else {
      const p1 = S === 1 ? loc.pos : locEnd.pos, p2 = S === 1 ? locEnd.pos : loc.pos;
      const span = x.start === x.end ? p1 : p1 + '_' + p2.replace(/^[cn]\./, '');
      if (x.type === 'del') res.hgvs = span + 'del';
      else if (x.type === 'delins') res.hgvs = span + 'delins' + altT;
      else if (x.type === 'dup' && x.dupStart != null) {
        const a = locate(t, x.dupStart).pos, b = locate(t, x.dupEnd).pos;
        const [u, w] = S === 1 ? [a, b] : [b, a];
        res.hgvs = (x.dupStart === x.dupEnd ? u : u + '_' + w.replace(/^[cn]\./, '')) + 'dup';
      } else if (x.type === 'ins' || x.type === 'dup') {
        const insT = S === 1 ? x.alt.slice(1) : rc(x.alt.slice(1));
        const a = locate(t, x.start).pos, b = locate(t, x.start + 1).pos;
        const [u, w] = S === 1 ? [a, b] : [b, a];
        res.hgvs = `${u}_${w.replace(/^[cn]\./, '')}ins${insT}`;
      }
    }
    const desc = { upstream: 'upstream of this isoform', downstream: 'downstream of this isoform', intron: 'intron', exon: 'exon', unknown: '' };
    res.where = loc.type === 'exon' ? `exon ${loc.exon + 1}/${t.exons.length}` : loc.type === 'intron' ? `intron ${loc.intron + 1}/${t.exons.length - 1}` :
      loc.type === 'upstream' ? `${fmt(loc.d)} bp upstream of the transcription start` : loc.type === 'downstream' ? `${fmt(loc.d)} bp downstream of the 3′ end` : desc[loc.type];
    if (loc.type === 'upstream' || loc.type === 'downstream') {
      res.consequence = loc.type + '_transcript_variant'; res.effect = 'not in this mRNA';
      res.summary = `Outside this isoform: ${res.where}. Its mRNA and protein sequences are unchanged.`;
      return res;
    }
    if (x.type === 'pos') { res.consequence = 'position only'; res.effect = loc.type === 'exon' ? 'exonic position' : 'intronic position'; return res; }
    // splice-site proximity (all variant types)
    if (loc.type === 'intron') {
      const d5 = loc.dUp, d3 = loc.dDown;
      if (d5 <= 2 || d3 <= 2) { res.consequence = d5 <= 2 ? 'splice_donor_variant' : 'splice_acceptor_variant'; }
      else if (d5 <= 8 || d3 <= 8) res.consequence = 'splice_region_variant';
      else if (d3 <= 17) res.consequence = 'splice_polypyrimidine_tract_variant';
      else res.consequence = 'intron_variant';
      res.effect = /donor|acceptor/.test(res.consequence) ? 'splice site disrupted' : 'intronic';
      res.summary = `Intronic (${loc.pos}); ${d5}/${d3} nt from the flanking exons.`;
      if (/donor|acceptor/.test(res.consequence)) {
        const skipK = d5 <= 2 ? loc.intron : loc.intron + 1; // donor: exon upstream; acceptor: exon downstream
        res.skip = skipScenario(v, t, skipK);
        res.summary = `Canonical splice ${d5 <= 2 ? 'donor' : 'acceptor'} (${loc.pos}). Likely outcomes: skipping of exon ${skipK + 1} (${res.skip.text}) or intron retention.`;
      }
      return res;
    }
    // exonic
    const k = loc.exon, i = loc.i, exLen = t.exons[k][1] - t.exons[k][0] + 1, offIn = i - G.cum[k];
    if ((k > 0 && offIn < 3) || (k < t.exons.length - 1 && exLen - 1 - offIn < 3)) res.splice_region = true;
    if (G.cds0 < 0) { res.consequence = 'non_coding_exon_variant'; res.effect = 'non-coding exon'; res.summary = `Non-coding exon ${k + 1}.`; return res; }
    const m = mrna(v, t);
    if (x.type === 'sub') {
      if (i < G.cds0) { res.consequence = '5_prime_UTR_variant'; res.effect = "5′ UTR"; res.summary = `5′ UTR of exon ${k + 1}.`;
        const mm = m.slice(0, i) + altT + m.slice(i + 1); const ctx = mm.slice(Math.max(0, i - 2), i + 3);
        if (/ATG/.test(ctx) && !/ATG/.test(m.slice(Math.max(0, i - 2), i + 3))) res.summary += ' Creates an upstream ATG.';
        return res; }
      if (i > G.cds1) { res.consequence = '3_prime_UTR_variant'; res.effect = "3′ UTR"; res.summary = `3′ UTR of exon ${k + 1}.`; return res; }
      const ci = i - G.cds0, q = Math.floor(ci / 3), p = ci % 3;
      const codon = m.substr(G.cds0 + 3 * q, 3), mut = codon.slice(0, p) + altT + codon.slice(p + 1);
      const a1 = CODON[codon] || 'X', a2 = CODON[mut] || 'X';
      Object.assign(res, { codon: codon + '>' + mut, codon_pos: p + 1, aa_ref: a1, aa_alt: a2, prot_pos: q + 1, cds_pos: ci + 1 });
      const protLen = (G.cds1 - G.cds0 + 1) / 3 - 1; res.prot_len = protLen;
      if (q === 0 && a1 === 'M' && a2 !== 'M') { res.consequence = 'start_lost'; res.hgvs_p = 'p.Met1?'; res.effect = 'start codon lost'; res.summary = 'Start codon lost: translation would begin at a downstream ATG, if at all.'; }
      else if (a1 === '*' && a2 !== '*') { const ext = translateFrom(m.slice(0, i) + altT + m.slice(i + 1), G.cds0 + 3 * q);
        res.consequence = 'stop_lost'; res.hgvs_p = `p.Ter${q + 1}${AA3[a2]}ext${ext.endsWith('*') ? 'Ter' + (ext.length - 1) : 'Ter?'}`; res.effect = 'stop lost'; res.summary = `Stop codon lost; the protein runs into the 3′ UTR${ext.endsWith('*') ? ` and gains ${ext.length - 1} residues` : ' without a new stop'}.`; }
      else if (a2 === '*') { res.consequence = 'stop_gained'; res.hgvs_p = `p.${AA3[a1]}${q + 1}Ter`; res.hgvs_p1 = `p.${a1}${q + 1}*`;
        res.ptc = { i: G.cds0 + 3 * q, prot_len: q }; }
      else if (a1 === a2) { res.consequence = 'synonymous_variant'; res.hgvs_p = `p.${AA3[a1]}${q + 1}=`; res.effect = 'synonymous'; res.summary = `Synonymous (${codon}>${mut}, ${AA3[a1]}${q + 1}).`; }
      else { res.consequence = 'missense_variant'; res.hgvs_p = `p.${AA3[a1]}${q + 1}${AA3[a2]}`; res.hgvs_p1 = `p.${a1}${q + 1}${a2}`; res.effect = 'missense'; res.summary = `Missense ${AA3[a1]}${q + 1}${AA3[a2]} (${codon}>${mut}).`; }
      if (res.ptc) finishPTC(v, t, res, t.exons.map(([s, e]) => e - s + 1));
      if (res.splice_region) res.summary = (res.summary || '') + ` Within 3 nt of an exon boundary: may also alter splicing.`;
      return res;
    }
    // indels within one exon
    const ka = exonOf(t, x.start), kb = exonOf(t, x.end);
    if (ka < 0 || kb < 0 || ka !== kb) { res.consequence = 'complex (crosses an exon boundary)'; res.effect = 'splice structure changed'; res.summary = 'The change crosses an exon-intron boundary; splicing is likely altered and is not modelled here.'; return res; }
    const ia = g2i(t, x.start), ib = g2i(t, x.end), i0 = Math.min(ia, ib), i1 = Math.max(ia, ib);
    let mm;
    if (x.type === 'ins' || x.type === 'dup') {
      const insT = S === 1 ? x.alt.slice(1) : rc(x.alt.slice(1));
      const at = S === 1 ? i0 + 1 : i0; // + strand: after anchor; - strand: before anchor (transcript orientation)
      mm = m.slice(0, at) + insT + m.slice(at);
    } else mm = m.slice(0, i0) + altT + m.slice(i1 + 1);
    const delta = mm.length - m.length;
    if (i1 < G.cds0 || i0 > G.cds1) { res.consequence = i1 < G.cds0 ? '5_prime_UTR_variant' : '3_prime_UTR_variant'; res.effect = i1 < G.cds0 ? '5′ UTR' : '3′ UTR'; res.summary = `UTR ${x.type} in exon ${k + 1}.`; return res; }
    const refP = translateFrom(m, G.cds0), newP = translateFrom(mm, G.cds0);
    let f = 0; while (f < refP.length && f < newP.length && refP[f] === newP[f]) f++;
    const lens = t.exons.map(([s, e]) => e - s + 1); lens[k] += delta;
    if (delta % 3 !== 0) {
      res.consequence = 'frameshift_variant';
      if (newP.endsWith('*')) { const stopIdx = newP.length - 1;
        res.hgvs_p = newP[f] === '*' ? `p.${AA3[refP[f]]}${f + 1}Ter` : `p.${AA3[refP[f]]}${f + 1}${AA3[newP[f]] || 'Ter'}fsTer${stopIdx - f + 1}`;
        res.ptc = { i: G.cds0 + 3 * stopIdx, prot_len: stopIdx }; }
      else { res.hgvs_p = `p.${AA3[refP[f]]}${f + 1}fs`; res.effect = 'frameshift, no stop before the polyA'; res.summary = 'Frameshift with no downstream stop codon (non-stop mRNA).'; }
      Object.assign(res, { prot_pos: f + 1, aa_ref: refP[f], prot_len: refP.length - 1 });
    } else if (delta === 0) {
      Object.assign(res, { prot_pos: f + 1, aa_ref: refP[f], prot_len: refP.length - 1 });
      let l = refP.length - 1; while (l > f && refP[l] === newP[l]) l--;
      const stopIdx = newP.indexOf('*');
      if (stopIdx >= 0 && stopIdx < refP.length - 1) { res.consequence = 'stop_gained'; res.hgvs_p = `p.${AA3[refP[stopIdx]]}${stopIdx + 1}Ter`; res.ptc = { i: G.cds0 + 3 * stopIdx, prot_len: stopIdx }; }
      else if (f >= refP.length) { res.consequence = 'synonymous_variant'; res.effect = 'synonymous'; res.summary = 'Multi-nucleotide change without an amino-acid change.'; }
      else { const span = l > f ? `${AA3[refP[f]]}${f + 1}_${AA3[refP[l]]}${l + 1}delins${[...newP.slice(f, l + 1)].map(a => AA3[a]).join('')}` : `${AA3[refP[f]]}${f + 1}${AA3[newP[f]]}`;
        res.consequence = 'missense_variant'; res.hgvs_p = 'p.' + span; res.effect = 'missense'; res.summary = `Multi-nucleotide substitution: ${res.hgvs_p}.`; }
    } else {
      const expected = (refP.length - 1) + delta / 3, stopIdx = newP.endsWith('*') ? newP.length - 1 : -1;
      Object.assign(res, { prot_pos: f + 1, aa_ref: refP[f], prot_len: refP.length - 1 });
      if (stopIdx >= 0 && stopIdx < expected) {
        res.consequence = 'stop_gained'; res.hgvs_p = `p.(${AA3[refP[f]] || ''}${f + 1}: in-frame change creates a stop at codon ${stopIdx + 1})`;
        res.ptc = { i: G.cds0 + 3 * stopIdx, prot_len: stopIdx };
      } else {
        res.consequence = delta < 0 ? 'inframe_deletion' : 'inframe_insertion';
        res.hgvs_p = `p.(${delta < 0 ? (-delta / 3) + ' aa deleted' : (delta / 3) + ' aa inserted'} at ${AA3[refP[f]] || ''}${f + 1})`;
        res.effect = 'in-frame'; res.summary = `In-frame ${delta < 0 ? 'deletion' : 'insertion'} of ${Math.abs(delta) / 3} residue(s) at ${AA3[refP[f]] || ''}${f + 1}; no premature stop.`;
      }
    }
    if (res.ptc) finishPTC(v, t, res, lens);
    return res;
  }
  function finishPTC(v, t, res, lens) {
    const G = geom(t);
    const r = nmd(lens, res.ptc.i, G.cds0); res.nmd = r;
    const doms = domainsFor(v, t); const L = res.ptc.prot_len;
    res.lost = doms.filter(d => d.end > L).map(d => (d.start > L ? '' : 'part of ') + d.label);
    res.kept_len = L;
    const lostTxt = summariseLost(doms, L);
    if (r.verdict === 'NMD') { res.effect = 'NMD: mRNA degraded'; res.summary = `Premature stop (${res.hgvs_p}). ${r.rule}: nonsense-mediated decay predicted, so little mRNA or protein is expected. Any escaping protein would stop at residue ${fmt(L)}${lostTxt ? ', losing ' + lostTxt : ''}.`; }
    else if (r.verdict === 'possible escape') { res.effect = 'may escape NMD'; res.summary = `Premature stop (${res.hgvs_p}); ${r.rule}. The mRNA may partly escape NMD and make a protein of ${fmt(L)} residues${lostTxt ? ' lacking ' + lostTxt : ''}.`; }
    else { res.effect = 'escapes NMD: truncated protein'; res.summary = `Premature stop (${res.hgvs_p}); ${r.rule}, so NMD is not expected. Predicted product: ${fmt(L)} residues${lostTxt ? ' lacking ' + lostTxt : ''}.`; }
    // a downstream in-frame Met (possible re-initiation)
    const prot = translate(mrna(v, t).slice(G.cds0, G.cds1 + 1)); const nm = prot.indexOf('M', L + 1);
    if (nm > 0) res.next_met = nm + 1;
  }
  function summariseLost(doms, L) {
    const lost = doms.filter(d => d.end > L && !['unique', 'extra'].includes(d.type)).map(d => d.type === 'dv23' ? Object.assign({}, d, { label: 'DV23' }) : d);
    if (!lost.length) return '';
    const srs = lost.filter(d => d.type === 'SR'), other = lost.filter(d => d.type !== 'SR');
    const parts = [];
    if (srs.length) { const nums = srs.map(d => +(d.label.match(/\d+/) || [0])[0]).filter(Boolean);
      parts.push(nums.length ? `spectrin repeats ${Math.min(...nums)}–${Math.max(...nums)}` : `${srs.length} spectrin repeats`); }
    for (const d of other) if (!parts.includes(d.label)) parts.push(d.label);
    return parts.join(', ');
  }
  function skipScenario(v, t, k) {
    const G = geom(t); const L = t.exons[k][1] - t.exons[k][0] + 1;
    if (G.cds0 < 0) return { text: 'non-coding isoform' };
    const cdsInExon = Math.max(0, Math.min(G.cum[k + 1] - 1, G.cds1) - Math.max(G.cum[k], G.cds0) + 1);
    if (!cdsInExon) return { text: 'exon outside the CDS' };
    if (k === 0 || k === t.exons.length - 1) return { text: 'terminal exon: skipping unlikely' };
    const m = mrna(v, t); const mm = m.slice(0, G.cum[k]) + m.slice(G.cum[k + 1]);
    const lens = t.exons.map(([s, e]) => e - s + 1); lens.splice(k, 1);
    const newP = translateFrom(mm, G.cds0);
    if (cdsInExon % 3 === 0 && newP.endsWith('*') && newP.length - 1 === (G.cds1 - G.cds0 + 1 - cdsInExon) / 3 - 1)
      return { text: `in frame, removes ${cdsInExon / 3} residues; no PTC` };
    const st = newP.length - 1; const r = nmd(lens, G.cds0 + 3 * st, G.cds0);
    return { text: `${cdsInExon % 3 ? 'frameshift' : 'in frame'}, new stop after ${fmt(st)} residues; ${r.verdict === 'NMD' ? 'NMD predicted' : r.verdict}`, nmd: r };
  }

  /* ------------------------------------------------------------------ orthologous positions */
  function orthologs(DATA, v, t, protPos) {
    const msa = DATA.alignments && DATA.alignments.giant_protein, vm = DATA.giant_msa_views || {};
    if (!msa || t.id.split('.')[0] !== String(v.giant_ref).split('.')[0]) return [];
    const rowOf = Object.keys(vm).find(k => vm[k] === v.key); if (!rowOf) return [];
    const row = msa.rows.find(r => r[0] === rowOf); if (!row) return [];
    let col = -1, n = 0; for (let c = 0; c < row[1].length; c++) if (row[1][c] !== '-') { n++; if (n === protPos) { col = c; break; } }
    if (col < 0) return [];
    const out = [];
    for (const [name, s] of msa.rows) {
      const key = vm[name]; if (!key) continue;
      const ov = (DATA.views || []).find(x => x.key === key); if (!ov) continue;
      const aa = s[col]; if (aa === '-') { out.push({ view: key, species: ov.species, symbol: ov.symbol, gap: true }); continue; }
      const pos = s.slice(0, col + 1).replace(/-/g, '').length;
      const gt = byId(ov, ov.giant_ref); const G = geom(gt); const ci = G.cds0 + 3 * (pos - 1);
      const g1 = i2g(gt, ci), g3 = i2g(gt, ci + 2); const codon = mrna(ov, gt).substr(ci, 3);
      out.push({ view: key, species: ov.species, symbol: ov.symbol, assembly: ov.assembly, chrom: ov.chrom, transcript: gt.id, aa, pos, codon,
        exon: exonOf(gt, g1) + 1, nExons: gt.exons.length, gStart: Math.min(g1, g3), gEnd: Math.max(g1, g3), strand: gt.strand,
        self: key === v.key });
    }
    return out;
  }

  /* ------------------------------------------------------------------ analyse (top level) */
  function analyse(DATA, v, input, opts) {
    opts = opts || {};
    const P = parse(input);
    if (P.error) return { error: P.error };
    const R = resolve(v, P, opts.tx);
    if (R.error) return { error: R.error, notes: R.notes || [] };
    const variants = R.variants.map(x => analyseOne(DATA, v, x));
    return { input, parsed: P, notes: R.notes, variants };
  }
  function analyseOne(DATA, v, x) {
    const g = { chrom: v.chrom, acc: v.chrom_acc, assembly: v.assembly, start: x.start, end: x.end, ref: x.ref, alt: x.alt, type: x.type, ref_ok: x.ref_ok, given_ref: x.given_ref };
    if (x.type === 'sub') g.hgvs_g = `${v.chrom_acc || 'chr' + v.chrom}:g.${x.start}${x.ref}>${x.alt}`;
    else if (x.type === 'del') g.hgvs_g = `${v.chrom_acc}:g.${x.start}${x.end > x.start ? '_' + x.end : ''}del`;
    else if (x.type === 'delins') g.hgvs_g = `${v.chrom_acc}:g.${x.start}${x.end > x.start ? '_' + x.end : ''}delins${x.alt}`;
    else if (x.type === 'ins' || x.type === 'dup') g.hgvs_g = `${v.chrom_acc}:g.${x.start}_${x.start + 1}ins${x.alt.slice(1)}`;
    g.vcf = x.type === 'sub' ? `${v.chrom}-${x.start}-${x.ref}-${x.alt}` : null;
    if (v.grch37_offset) g.grch37 = { start: x.start + v.grch37_offset, end: x.end + v.grch37_offset };
    g.context = v.seq.slice(Math.max(0, x.start - v.region.start - 12), x.start - v.region.start).toLowerCase() + '[' + (x.ref || '') + (x.alt != null ? '>' + x.alt : '') + ']' +
      v.seq.slice(x.end - v.region.start + 1, x.end - v.region.start + 13).toLowerCase();
    const per = v.transcripts.map(t => consequence(v, t, x));
    const giant = per.find(r => r.id === v.giant_ref);
    const a2 = per.filter(r => r.category === 'alpha2');
    const focus = per.find(r => r.id === x.tx) || giant;
    let orth = [];
    if (giant && giant.prot_pos && (giant.consequence === 'stop_gained' || /missense|synonymous|frameshift|inframe|stop_lost|start_lost/.test(giant.consequence || '')))
      orth = orthologs(DATA, v, byId(v, v.giant_ref), giant.prot_pos);
    // relation to the alpha2 promoter / first exon
    let a2rel = null;
    const a2t = v.transcripts.find(t => t.category === 'alpha2');
    if (a2t) { const G = geom(a2t); const d = (x.start - G.tss) * -a2t.strand; // >0: upstream of TSS in transcript orientation
      a2rel = { tx: a2t.id, upstream_of_tss: d > 0 ? d : null, inside: d <= 0 && exonOf(a2t, x.start) >= 0 };
    }
    return { genomic: g, per, giant, alpha2: a2, focus, orthologs: orth, alpha2_relation: a2rel, from_protein: x.from_protein || null };
  }

  /* ------------------------------------------------------------------ figure */
  const C = { ink: '#1b2430', ink2: '#55606e', ink3: '#8a95a3', line: '#cfd6de', bg: '#ffffff', panel: '#f5f7f9',
    crimson: '#b0124c', crimsonSoft: '#f4cbd9', teal: '#0a7b74', tealSoft: '#bfe2dd', slate: '#607083', slateSoft: '#d7dee6',
    red: '#c62828', amber: '#a15c00', green: '#1d7a4c', blue: '#2a5bd7', purple: '#7b3fa0', kash: '#222831', srA: '#d9dee5', srB: '#c3cad3',
    dv23: '#e0a33a', extra: '#9fb0c2' };
  const FONT = 'Helvetica, Arial, sans-serif', MONO = 'Menlo, Consolas, monospace';
  function T(x, y, s, o) { o = o || {}; return `<text x="${x.toFixed(1)}" y="${y.toFixed(1)}" font-family="${o.mono ? MONO : FONT}" font-size="${o.size || 12}"${o.weight ? ` font-weight="${o.weight}"` : ''} fill="${o.fill || C.ink}"${o.anchor ? ` text-anchor="${o.anchor}"` : ''}${o.halo ? ' stroke="#ffffff" stroke-width="3" paint-order="stroke"' : ''}>${esc(s)}</text>`; }
  function shortRule(r) {
    if (r.nmd) return r.nmd.verdict === 'NMD' ? `PTC ${fmt(r.nmd.distToLastEJ)} nt before the last exon junction · exon ${r.nmd.exon}/${r.nmd.nExons}` : r.nmd.rule;
    return r.where || '';
  }
  function verdictStyle(r) {
    if (!r || r.effect === 'not in this mRNA') return { fill: C.green, text: 'variant not in this mRNA: unaffected' };
    if (r.nmd && r.nmd.verdict === 'NMD') return { fill: C.red, text: 'PTC → NMD: mRNA degraded' };
    if (r.nmd && r.nmd.verdict === 'possible escape') return { fill: C.amber, text: 'PTC: may escape NMD' };
    if (r.nmd) return { fill: C.amber, text: 'PTC escapes NMD: truncated protein' };
    if (/intron|UTR|synonymous/.test(r.consequence || '')) return { fill: C.ink2, text: r.effect || r.consequence };
    return { fill: C.amber, text: r.effect || r.consequence };
  }
  function pickRows(v, A) {
    const rows = [];
    const add = r => { if (r && !rows.includes(r) && rows.length < 5) rows.push(r); };
    add(A.giant); (A.alpha2 || []).slice(0, 1).forEach(add);
    const rank = r => (r.ptc ? 0 : r.consequence && /missense|frameshift|splice_donor|splice_acceptor|inframe/.test(r.consequence) ? 1 : 2) * 10 + (/^NM_/.test(r.id) ? 0 : /^XM_/.test(r.id) ? 1 : 2);
    A.per.filter(r => r !== A.giant && r.category !== 'alpha2' && r.effect !== 'not in this mRNA' && r.loc.type === 'exon').sort((a, b) => rank(a) - rank(b)).forEach(add);
    return rows;
  }
  function figureSVG(DATA, v, A, opts) {
    opts = opts || {};
    const W = 1100, L = 210, R = W - 28, PW = R - L;
    const g = A.genomic, giantT = byId(v, v.giant_ref), rows = pickRows(v, A).map(r => ({ r, t: byId(v, r.id) }));
    const varLabel = A.giant && A.giant.hgvs ? `${v.giant_ref}:${A.giant.hgvs}${A.giant.hgvs_p ? ' (' + A.giant.hgvs_p + ')' : ''}` : (g.hgvs_g || '');
    let y = 0, s = '';
    // header
    s += T(24, 32, `${v.symbol} ${varLabel}`, { size: 19, weight: 700 });
    const refTxt = g.ref_ok === false ? 'REFERENCE MISMATCH' : 'reference allele verified';
    s += T(24, 54, `${v.species} ${v.assembly} chr${v.chrom}:${fmt(g.start)}${g.end > g.start ? '-' + fmt(g.end) : ''} ${g.ref || ''}${g.alt != null ? '>' + (g.alt || '-') : ''} (+ strand)` +
      (g.grch37 ? ` · GRCh37 chr${v.chrom}:${fmt(g.grch37.start)}` : '') + ` · ${refTxt}`, { size: 12.5, fill: g.ref_ok === false ? C.red : C.ink2 });
    y = 78;
    // ---------- panel A: whole gene in transcript orientation
    const orient = giantT.strand, G0 = geom(giantT);
    const gA = Math.min(G0.lo, g.start - 1000), gB = Math.max(G0.hi, g.end + 1000);
    const sxA = gpos => orient === 1 ? L + (gpos - gA) / (gB - gA) * PW : L + (gB - gpos) / (gB - gA) * PW;
    s += T(24, y + 4, 'A', { size: 15, weight: 700 }) + T(44, y + 4, `Gene (${v.symbol}, 5′→3′, ${fmt(gB - gA)} bp)`, { size: 13, weight: 700 });
    y += 22;
    s += `<line x1="${L}" y1="${y + 14}" x2="${R}" y2="${y + 14}" stroke="${C.ink3}" stroke-width="1.2"/>`;
    for (const [a, b] of giantT.exons) { const x1 = sxA(a), x2 = sxA(b); s += `<rect x="${Math.min(x1, x2).toFixed(1)}" y="${y + 8}" width="${Math.max(1, Math.abs(x2 - x1)).toFixed(1)}" height="12" fill="${C.teal}"/>`; }
    s += T(24, y + 19, `${v.giant_ref} (giant)`, { size: 11.5, fill: C.teal, weight: 600 });
    const a2t = rows.map(o => o.t).find(t => t && t.category === 'alpha2') || v.transcripts.find(t => t.category === 'alpha2');
    if (a2t) { const Ga = geom(a2t); const x1 = sxA(Ga.lo), x2 = sxA(Ga.hi);
      s += `<path d="M${Math.min(x1, x2).toFixed(1)} ${y + 3} v-5 h${Math.abs(x2 - x1).toFixed(1)} v5" fill="none" stroke="${C.crimson}" stroke-width="1.5"/>`;
      const lab = `${v.symbol} α2 (${a2t.predicted ? 'predicted' : a2t.id})`, lw = lab.length * 6.3, cx = (x1 + x2) / 2;
      s += T(Math.min(Math.max(cx - lw / 2, L), R - lw), y - 6, lab, { size: 11, fill: C.crimson, weight: 600 }); }
    const vx = sxA(g.start);
    s += `<line x1="${vx.toFixed(1)}" y1="${y - 2}" x2="${vx.toFixed(1)}" y2="${y + 30}" stroke="${C.red}" stroke-width="2"/>` +
      `<path d="M${(vx - 6).toFixed(1)} ${y + 38} l6 -8 l6 8 z" fill="${C.red}"/>` + T(vx, y + 52, A.giant ? `${A.giant.hgvs} · ${A.giant.where}` : g.hgvs_g, { size: 11.5, fill: C.red, anchor: vx > R - 160 ? 'end' : vx < L + 160 ? 'start' : 'middle', weight: 600 });
    s += T(L, y + 52, '5′', { size: 11, fill: C.ink3 }) + T(R, y + 52, '3′', { size: 11, fill: C.ink3, anchor: 'end' });
    y += 76;
    // ---------- panel B: zoom on variant and alpha2 start
    let zA, zB;
    const near = [g.start];
    if (a2t) { near.push(a2t.exons[0][0], a2t.exons[0][1]); if (a2t.exons[1]) near.push(a2t.exons[1][0], a2t.exons[1][1]); }
    zA = Math.min(...near); zB = Math.max(...near);
    if (zB - zA > 60000 || zB - zA < 2000) { zA = g.start - 4000; zB = g.end + 4000; }
    const pad = Math.max(600, (zB - zA) * 0.08); zA -= pad; zB += pad;
    const sxB = gpos => orient === 1 ? L + (gpos - zA) / (zB - zA) * PW : L + (zB - gpos) / (zB - zA) * PW;
    s += T(24, y + 4, 'B', { size: 15, weight: 700 }) + T(44, y + 4, `Variant and neighbouring exons (${fmt(Math.round(zB - zA))} bp window)`, { size: 13, weight: 700 });
    y += 30;
    const yTop = y, exLabels = [];
    rows.forEach(({ r, t }, j) => {
      const cy = y + 12, col = t.category === 'alpha2' ? C.crimson : t.category === 'giant' ? C.teal : C.slate;
      s += T(24, cy + 4, (t.predicted ? t.name : t.id).slice(0, 30), { size: 11.5, fill: col, weight: 600, mono: true });
      const Gt = geom(t); const lo = Math.max(Gt.lo, zA), hi = Math.min(Gt.hi, zB);
      if (hi > lo) s += `<line x1="${sxB(lo).toFixed(1)}" y1="${cy}" x2="${sxB(hi).toFixed(1)}" y2="${cy}" stroke="${col}" stroke-width="1"/>`;
      t.exons.forEach(([a, b], k) => {
        if (b < zA || a > zB) return;
        const x1 = sxB(Math.max(a, zA)), x2 = sxB(Math.min(b, zB)); const xl = Math.min(x1, x2), w = Math.max(2, Math.abs(x2 - x1));
        const unique = t.category === 'alpha2' && k === 0;
        const coding = t.cds && b >= t.cds[0] && a <= t.cds[1];
        const h = coding ? 14 : 8;
        s += `<rect x="${xl.toFixed(1)}" y="${(cy - h / 2).toFixed(1)}" width="${w.toFixed(1)}" height="${h}" fill="${unique ? C.crimson : t.category === 'alpha2' ? C.crimsonSoft : col}" stroke="${t.category === 'alpha2' ? C.crimson : 'none'}"${t.predicted && unique ? ' stroke-dasharray="3 2"' : ''}/>`;
        exLabels.push(T(xl + w / 2, cy - 10, unique ? '1α2' : 'ex ' + (k + 1), { size: 10, fill: unique ? C.crimson : C.ink2, anchor: 'middle', halo: true }));
      });
      y += 36;
    });
    const vxB = sxB(g.start);
    s += `<line x1="${vxB.toFixed(1)}" y1="${yTop - 8}" x2="${vxB.toFixed(1)}" y2="${y - 6}" stroke="${C.red}" stroke-width="1.5" stroke-dasharray="4 3"/>` + exLabels.join('');
    s += T(vxB, yTop - 12, 'variant', { size: 11, fill: C.red, anchor: 'middle', weight: 600 });
    if (A.alpha2_relation && A.alpha2_relation.upstream_of_tss && a2t) {
      const xT = sxB(geom(a2t).tss), yb = y + 2;
      s += `<path d="M${vxB.toFixed(1)} ${yb} v6 H${xT.toFixed(1)} v-6" fill="none" stroke="${C.crimson}" stroke-width="1.2"/>`;
      s += T((vxB + xT) / 2, yb + 20, `${fmt(A.alpha2_relation.upstream_of_tss)} bp from the variant to the α2 transcription start (exon 1α2)`, { size: 11.5, fill: C.crimson, anchor: 'middle' });
      y += 30;
    }
    y += 16;
    // ---------- panel C: mRNA fate
    s += T(24, y + 4, 'C', { size: 15, weight: 700 }) + T(44, y + 4, 'mRNA fate (each mRNA scaled to its own length)', { size: 13, weight: 700 });
    y += 26;
    const BW = 450, bx = L;
    rows.forEach(({ r, t }) => {
      const Gt = geom(t), sc = BW / Gt.len, cy = y + 12;
      const col = t.category === 'alpha2' ? C.crimson : t.category === 'giant' ? C.teal : C.slate;
      s += T(24, cy + 4, (t.predicted ? t.name : t.id).slice(0, 30), { size: 11.5, fill: col, weight: 600, mono: true });
      t.exons.forEach(([a, b], k) => {
        const x1 = bx + Gt.cum[k] * sc, w = Math.max(0.6, (Gt.cum[k + 1] - Gt.cum[k]) * sc);
        const cdsA = Math.max(Gt.cum[k], Gt.cds0), cdsB = Math.min(Gt.cum[k + 1] - 1, Gt.cds1);
        s += `<rect x="${x1.toFixed(2)}" y="${cy - 3}" width="${w.toFixed(2)}" height="6" fill="${k % 2 ? C.slateSoft : C.srB}"/>`;
        if (Gt.cds0 >= 0 && cdsA <= cdsB) s += `<rect x="${(bx + cdsA * sc).toFixed(2)}" y="${cy - 6}" width="${Math.max(0.6, (cdsB - cdsA + 1) * sc).toFixed(2)}" height="12" fill="${t.category === 'alpha2' ? (k === 0 ? C.crimson : C.crimsonSoft) : k % 2 ? C.slateSoft : C.srB}" stroke="${t.category === 'alpha2' && k ? C.crimson : 'none'}" stroke-width="0.4"/>`;
      });
      const lastEJ = bx + Gt.cum[t.exons.length - 1] * sc;
      s += `<line x1="${lastEJ.toFixed(1)}" y1="${cy - 11}" x2="${lastEJ.toFixed(1)}" y2="${cy + 11}" stroke="${C.ink}" stroke-width="1"/>`;
      if (r.ptc) { const px = bx + r.ptc.i * sc;
        s += `<polygon points="${[0, 1, 2, 3, 4, 5, 6, 7].map(q => { const a = Math.PI / 8 + q * Math.PI / 4; return (px + 7 * Math.cos(a)).toFixed(1) + ',' + (cy - 18 + 7 * Math.sin(a)).toFixed(1); }).join(' ')}" fill="${C.red}"/>` +
          `<line x1="${px.toFixed(1)}" y1="${cy - 11}" x2="${px.toFixed(1)}" y2="${cy + 8}" stroke="${C.red}" stroke-width="2"/>`; }
      else if (r.loc.type === 'exon' || r.loc.type === 'intron') { const ii = r.loc.type === 'exon' ? r.loc.i : Gt.cum[r.loc.intron + 1]; const px = bx + ii * sc;
        s += `<path d="M${(px - 5).toFixed(1)} ${cy - 16} l5 7 l5 -7 z" fill="${C.amber}"/>`; }
      const vs = verdictStyle(r);
      s += `<rect x="${bx + BW + 16}" y="${cy - 10}" width="9" height="20" fill="${vs.fill}"/>` + T(bx + BW + 32, cy - 1, vs.text, { size: 12, fill: vs.fill, weight: 700 }) +
        T(bx + BW + 32, cy + 13, shortRule(r), { size: 10.5, fill: C.ink2 });
      y += 44;
    });
    s += T(L, y - 6, 'vertical black tick = last exon-exon junction; red stop sign = premature stop codon; thick = coding', { size: 10.5, fill: C.ink3 });
    y += 18;
    // ---------- panel D: protein
    s += T(24, y + 4, 'D', { size: 15, weight: 700 }) + T(44, y + 4, 'Protein domains (each protein scaled to its own length; UniProt Q8NF91 domains projected through shared exons)', { size: 13, weight: 700 });
    y += 26;
    rows.forEach(({ r, t }) => {
      const doms = domainsFor(v, t), plen = (t.protein && t.protein.length) || r.prot_len || ((geom(t).cds1 - geom(t).cds0 + 1) / 3 - 1);
      if (!plen || plen < 1) return;
      const sc = BW / plen, cy = y + 12, col = t.category === 'alpha2' ? C.crimson : t.category === 'giant' ? C.teal : C.slate;
      s += T(24, cy - 1, (t.predicted ? t.name : t.protein_id || t.id).slice(0, 27), { size: 11, fill: col, weight: 600, mono: true }) +
        T(24, cy + 12, `${fmt(plen)} aa`, { size: 10.5, fill: C.ink2 });
      s += `<rect x="${bx}" y="${cy - 4}" width="${BW}" height="8" fill="#eef1f4" stroke="${C.line}" stroke-width="0.5"/>`;
      let sr = 0;
      for (const d of doms) {
        if (d.type === 'motif') continue;
        const x1 = bx + (d.start - 1) * sc, w = Math.max(0.8, (d.end - d.start + 1) * sc);
        const f = d.type === 'CH' ? C.blue : d.type === 'SR' ? (sr++ % 2 ? C.srB : C.srA) : d.type === 'KASH' ? C.kash : d.type === 'TM' ? '#000' : d.type === 'unique' ? C.crimson : d.type === 'dv23' ? C.dv23 : d.type === 'extra' ? C.extra : C.slateSoft;
        const h = d.type === 'SR' ? 12 : 16;
        s += `<rect x="${x1.toFixed(2)}" y="${(cy - h / 2).toFixed(1)}" width="${w.toFixed(2)}" height="${h}" fill="${f}"/>`;
      }
      for (const d of doms) if (d.type === 'motif') { const x1 = bx + (d.start - 1) * sc; s += `<rect x="${(x1 - 1).toFixed(1)}" y="${cy - 10}" width="2.5" height="20" fill="${C.purple}"/>`; }
      if (r.ptc) { const xL = bx + r.ptc.prot_len * sc;
        s += `<rect x="${xL.toFixed(1)}" y="${cy - 11}" width="${(bx + BW - xL).toFixed(1)}" height="22" fill="#ffffff" fill-opacity="0.72"/>` +
          `<line x1="${xL.toFixed(1)}" y1="${cy - 14}" x2="${xL.toFixed(1)}" y2="${cy + 14}" stroke="${C.red}" stroke-width="2"/>`;
        s += T(bx + BW + 16, cy - 1, `truncated after residue ${fmt(r.ptc.prot_len)}`, { size: 12, fill: C.red, weight: 700 }) +
          T(bx + BW + 16, cy + 13, r.lost && r.lost.length ? 'lost: ' + summariseLost(doms, r.ptc.prot_len) : '', { size: 10.5, fill: C.ink2 }); }
      else if (r.prot_pos && r.aa_alt) { const xL = bx + (r.prot_pos - 1) * sc; s += `<path d="M${(xL - 5).toFixed(1)} ${cy - 16} l5 7 l5 -7 z" fill="${C.amber}"/>` + T(bx + BW + 16, cy + 4, `${r.hgvs_p}`, { size: 12, fill: C.amber, weight: 700 }); }
      else s += T(bx + BW + 16, cy + 4, r.effect === 'not in this mRNA' ? 'full-length protein, sequence unchanged' : (r.effect || ''), { size: 12, fill: r.effect === 'not in this mRNA' ? C.green : C.ink2, weight: 600 });
      y += 40;
    });
    // legend
    y += 4;
    const lg = [[C.blue, 'CH (actin-binding)'], [C.srA, 'spectrin repeat'], [C.purple, 'LEWD (kinesin-binding)'], [C.kash, 'KASH + TM'], [C.crimson, 'α2-specific (exon 1α2)'], [C.dv23, 'DV23 exon'], [C.extra, 'other isoform-specific'], [C.teal, 'giant exons']];
    let lx = 24; for (const [f, txt] of lg) { s += `<rect x="${lx}" y="${y - 9}" width="14" height="10" fill="${f}"/>` + T(lx + 19, y, txt, { size: 11, fill: C.ink2 }); lx += 28 + txt.length * 5.9; }
    y += 20;
    s += T(24, y, 'NMD rule: a PTC >50–55 nt upstream of the last exon-exon junction triggers decay; exceptions: last exon, <150 nt after the start codon, exons >407 nt (Lindeboom et al. 2016). Predictions, not measurements.', { size: 10.5, fill: C.ink3 });
    y += 16;
    s += T(24, y, `Data: ${opts.source || 'Ensembl / NCBI RefSeq / UniProt'} · SYNE1 isoform browser`, { size: 10.5, fill: C.ink3 });
    y += 14;
    return `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${y}" viewBox="0 0 ${W} ${y}"><rect width="${W}" height="${y}" fill="${C.bg}"/>${s}</svg>`;
  }

  return { parse, resolve, analyse, analyseOne, consequence, nmd, domainsFor, orthologs, figureSVG, locate, translate, rc, AA3, geom, mrna, byId };
}));
