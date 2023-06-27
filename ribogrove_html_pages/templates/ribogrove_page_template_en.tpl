<h1 style="display:none">RiboGrove</h1>
<img src="/cager/img/RiboGrove_logo.png"><br><br>
<h2>Contents</h2>
<ul>
<li><a href="#overview">What is RiboGrove</a>
<ul class="ribogrove-nested-list">
<li><a href="#other-databases">RiboGrove and other 16S rRNA databases</a></li>
<li><a href="#genome-categories">Genome categories</a></li>
</ul>
</li>
<li><a href="#downloads">Downloads</a>
<ul class="ribogrove-nested-list">
  <li><a href="#current-release">Latest RiboGrove release — {{ ribogrove_release_number }}</li>
  <li><a href="#release-archive">RiboGrove release archive</li>
  <li><a href="#release-notes">Release notes</li>
</ul>
</li>
<li><a href="#summary">Statistical summary</a>
<ul class="ribogrove-nested-list">
<li><a href="#db-size">RiboGrove size</a></li>
<li><a href="#gene-lens">16S rRNA gene lengths</a></li>
<li><a href="#copy-number">16S rRNA gene copy number</a></li>
<li><a href="#top-longest">Top-10 longest 16S rRNA genes</a></li>
<li><a href="#top-shortest">Top-10 shortest 16S rRNA genes</a></li>
<li><a href="#top-copy-num">Top-10 genomes with the largest 16S rRNA copy numbers</a></li>
<li><a href="#top-var">Top-10 genomes with the highest intragenomic variability of 16S rRNA genes</a></li>
<li><a href="#primers-coverages">Coverage of primer pairs for different V-regions of bacterial 16S rRNA genes</a></li>
</ul>
</li>
<li><a href="#searching-data">Searching data in RiboGrove</a>
<ul class="ribogrove-nested-list">
<li><a href="#header-format">Header format</a></li>
<li><a href="#select-by-header">Search sequences by header</a></li>
<li><a href="#select-by-len">Search sequences by length</a></li>
<li><a href="#select-head">Selecting header data</a></li>
</ul>
</li>
<li><a href="#contacts">Contacts</a></li>
<li><a href="#citing-ribogrove">Citing RiboGrove</a></li>
<li><a href="#faq">Questions people ask about RiboGrove</a></li>
</ul>
<hr>
<div id="overview" class="pad-anchor"></div>
<h2>What is RiboGrove</h2>
<p>RiboGrove is a database of 16S rRNA gene sequences of bacteria and archaea.</p>
<p>RiboGrove is based on the <a href="https://ncbi.nlm.nih.gov/refseq/">RefSeq</a> database. It contains only full-length sequences of 16S rRNA genes, and the sequences are derived from completely assembled prokaryotic genomes deposited in RefSeq. Hence we posit high reliability of RiboGrove sequences.</p>
<div id="other-databases" class="pad-anchor"></div>
<h3>RiboGrove and other 16S rRNA databases</h3>
<p>Here is a summary showing what is the (qualitative) difference between RiboGrove and similar rRNA sequence databases, namely <a href="https://rrndb.umms.med.umich.edu/"><i>rrn</i>DB</a>, <a href="https://www.arb-silva.de/">Silva</a>, <a href="http://rdp.cme.msu.edu/">RDP</a>, and <a href="https://greengenes.secondgenome.com/">Greengenes</a>. Briefly, RiboGrove is inferior in sequence amount and diversity, but superior in sequence reliability.</p>
<table class="sum-table">
<tbody>
  <tr class="sumtab-row"><th></th><th class="alcentr">RiboGrove</th><th class="alcentr"><i>rrn</i>DB</th><th class="alcentr">Silva</th><th class="alcentr">RDP</th><th class="alcentr">Greengenes</th></tr>
  <tr class="sumtab-row"><td class="feature-col">Represented organisms</td><td class="alcentr">Bacteria<br>Archaea</td><td class="alcentr">Bacteria<br>Archaea</td><td class="alcentr">Bacteria<br>Archaea<br>Eukaryotes</td><td class="alcentr">Bacteria<br>Archaea<br>Eukaryotes</td><td class="alcentr">Bacteria<br>Archaea</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Represented ribosome subunits</td><td class="alcentr">Small</td><td class="alcentr">Small</td><td class="alcentr">Large<br>Small</td><td class="alcentr">Large<br>Small</td><td class="alcentr">Small</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Contains sequences from assembled genomes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Contains amplicon sequences</td><td class="alcentr">No</td><td class="alcentr">No</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Contains partial gene sequences</td><td class="alcentr">No</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Discriminates genome categories</td><td class="alcentr">Yes</td><td class="alcentr">No</td><td class="alcentr"><i>Not<br>applicable</i></td><td class="alcentr"><i>Not<br>applicable</i></td><td class="alcentr"><i>Not<br>applicable</i></td></tr>
</tbody>
</table>
<div id="genome-categories" class="pad-anchor"></div>
<h3>Genome categories</h3>
<p>All genomes used for RiboGrove construction were divided into three categories according to their expected reliability:</p>
<ol style="list-style:none;">
  <li><strong>Category 1</strong> (the highest reliability). Genomes showing no signs of a low-quality assembly and sequenced either with PacBio technology or with a combination “Oxford&nbsp;Nanopore&nbsp;+&nbsp;Illumina”.</li>
  <li><strong>Category 2.</strong> Genomes showing no signs of a low-quality assembly and sequenced with any other technology (or the technology is not specified).</li>
  <li><strong>Category 3</strong> (the lowest reliability). Genomes showing at least one sign of a low-quality assembly.</li>
</ol>
<p>Signs of a low-quality assembly are the following:</p>
<ul>
  <li>The genome contains degenerate base(s) in 16S rRNA gene sequences.</li>
  <li>The assembly includes at least one RefSeq record whose title contains the phrase “map unlocalized” and this record contains a 16S rRNA gene or a part of it.</li>
</ul>
<br>
<p>The software used for the RiboGrove construction can be found in the following GitHub repository: <a href="https://github.com/masikol/ribogrove-tools">ribogrove-tools</a>.</p>
<hr>
<div id="downloads" class="pad-anchor"></div>
<h2>Downloads</h2>
<div id="current-release" class="pad-anchor"></div>
<h3>Latest RiboGrove release — {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
<p>The release is based on RefSeq release {{ refseq_release }}.</p>
<ul>
<li>A fasta file of full-length 16S gene sequences. <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Download (gzipped fasta file, {{ final_fasta_fsize_fmt }} MB)</a></li>
<li>Metadata <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/metadata_ribogrove_{{ ribogrove_release_number }}.zip">Download (zip archive {{ metadata_fsize_fmt }} MB)</a><br>
  <details>
    <summary><small>What information exactly does the metadata contain?</small></summary>
    <small>
      <p>The metadata consists of the following files:</p>
      <ol>
      <li>
        <span class="samp-highl">discarded_sequences.fasta.gz</span><br>
        This is a fasta file of sequences, which were present in source RefSeq genomes and were annotated a 16S rRNA genes but which have been discarded according to their incompleteness, internal repeats etc. and thus haven’t been included into RiboGrove.
      </li>
      <li>
        <span class="samp-highl">source_RefSeq_genomes.tsv</span><br>
        This is a TSV files, which contains information about what genomes were used for the RiboGrove construction.
      </li>
      <li>
        <span class="samp-highl">gene_seqs_statistics.tsv</span>, <span class="samp-highl">discarded_gene_seqs_statistics.tsv</span><br>
        These are TSV files, which contain nucleotide conposition, size, genomic and taxonomic affiliation of the gene sequences. The first file describes final RiboGrove gene sequences, the second file describes discarded sequences.
      </li>
      <li>
        <span class="samp-highl">categories.tsv</span><br>
        This is a TSV file, which contains information about what genome categories were assigned to each genome and why. Moreover, it contains information about what sequencing technology was used to sequence each genome.
      </li>
      <li>
        <span class="samp-highl">taxonomy.tsv</span><br>
        This is a TSV file, which contains taxonomic affiliation of each genome and gene.
      </li>
      <li>
        <span class="samp-highl">intragenic_repeats.tsv</span><br>
        This is a TSV file, which contains information about intragenomic repeats found in gene sequences using <a href="https://github.com/deprekate/RepeatFinder">RepeatFinder</a>.
      </li>
      <li>
        <span class="samp-highl">entropy_summary.tsv</span><br>
        This is a TSV file, which contains summary of instragenomic variability of the 16S rRNA genes. Intragenomic variability are calculated only for the category&nbsp;1 genomes having more than one 16S rRNA gene. Intragenomic variability is evaluated using Shannon entropy. We align gene sequences from each genome using <a href="https://drive5.com/muscle5/">MUSCLE</a>, and then we calculate Shannon entropy for each multiple alignment column (i.e. base).
      </li>
      <li>
        <span class="samp-highl">primer_pair_genomic_coverage.tsv</span><br>
        This is a TSV file which contains genomic coverage of primer pairs targeting different V-regions of 16S rRNA genes. For example, for <i>Enterobacteriaceae</i>, genomic coverage of a primer pair is the per cent of <i>Enterobacteriaceae</i> genomes which contain at least one 16S rRNA gene that can (theoretically) produce a PCR product using the primer pair.
      </li>
    </ol></small>
  </details>
</li>
</ul>
<p>The fasta file is compressed with gzip, and the metadata file is a zip archive. To uncompress them, Linux and Mac OS users may use gzip and zip programs, they should be built-in. For Windows users, the free and open-source (de)compression program <a href="https://www.7-zip.org/">7-Zip</a> is available.</p>
<div id="release-archive" class="pad-anchor"></div>
<h3>RiboGrove release archive</h3>
<p>You can find all releases in the <a href="/cager/en/content/60-ribogrove-archive">RiboGrove release archive</a>.</p>
<div id="release-notes" class="pad-anchor"></div>
<h3>Release notes</h3>
<p><i>No important differences from the previous release.</i></p>
<p>You can find notes to all RiboGrove releases on the <a href="/cager/en/content/66-ribogrove-release-notes">release notes page</a>.</p>
<hr>
<div id="summary" class="pad-anchor"></div>
<h2>Statistical summary</h2>
<div id="db-size" class="pad-anchor"></div>
<table class="sum-table"><caption>RiboGrove size
</caption>
<tbody>
<tr><th></th><th class="numcol">Bacteria</th><th class="numcol">Archaea</th><th class="numcol">Total</th></tr>
<tr class="sumtab-row">
<td>Number of gene sequences</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Number of unique gene sequences</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Number of species</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Number of genomes</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Number of genomes of category 1</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Number of genomes of category 2</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Number of genomes of category 3</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Total'] }}</td>
</tr>
</tbody>
</table>
<div id="gene-lens" class="pad-anchor"></div>
<table class="sum-table"><caption>16S rRNA gene lengths
</caption>
<tbody>
<tr><th></th><th class="numcol">Bacteria</th><th class="numcol">Archaea</th></tr>
<tr class="sumtab-row">
<td>Minimum (bp)</td>
<td class="numcol">{{ ribogrove_len_dict['min']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['min']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>25th percentile (bp) <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['25perc']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['25perc']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Median (bp) <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['median']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['median']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>75th percentile (bp) <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['75perc']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['75perc']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Average (bp) <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['mean']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['mean']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Mode (bp) <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['modes']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['modes']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Maximum (bp)</td>
<td class="numcol">{{ ribogrove_len_dict['max']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['max']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Standard deviation (bp) <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['std']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['std']['Archaea'] }}</td>
</tr>
</tbody>
</table>
<p><sup>*</sup> Metrics marked with an asterisk were calculated with preliminary normalization, i.e. median within-species gene length was used for the summary.</p>
<div id="copy-number" class="pad-anchor"></div>
<table class="sum-table">
<caption>16S rRNA gene copy number</caption>
<tbody>
<tr><th class="numcol">Copy number <sup>*</sup></th><th class="alncenter" colspan="2">Bacteria</th><th class="alncenter" colspan="2">Archaea</th></tr>
<tr><td></td><td>Number of species</td><td>Per cent of species (%)</td><td>Number of species</td><td>Per cent of species (%)</td></tr>

{% for _, row in ribogrove_copy_number_df.iterrows() %}
<tr class="sumtab-row">
<td class="numcol">{{ row['copy_number'] }}</td>
<td class="numcol">{{ row['number_of_species_bacteria'] }}</td>
<td class="numcol">{{ row['percent_of_species_bacteria'] }}</td>
<td class="numcol">{{ row['number_of_species_archaea'] }}</td>
<td class="numcol">{{ row['percent_of_species_archaea'] }}</td>
</tr>
{% endfor %}

</tbody>
</table>
<p><sup>*</sup> These are median within-species copy numbers.</p>

<div id="top-longest" class="pad-anchor"></div>
<table class="sum-table"><caption>Top-10 longest 16S rRNA genes
</caption>
<tbody>
<tr><th class="alnleft">Organism</th><th class="numcol">Gene length (bp)</th><th class="alnleft">RiboGrove Sequence ID(s)</th><th class="alnleft">Assembly<br>accession</th></tr>
<tr>
<td colspan="4" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_longest_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_longest_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<div id="top-shortest" class="pad-anchor"></div>
<table class="sum-table"><caption>Top-10 shortest 16S rRNA genes
</caption>
<tbody>
<tr><th class="alnleft">Organism</th><th class="numcol">Gene length (bp)</th><th class="alnleft">RiboGrove Sequence ID(s)</th><th class="alnleft">Assembly<br>accession</th></tr>
<tr>
<td colspan="4" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_shortest_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_shortest_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<div id="top-copy-num" class="pad-anchor"></div>
<table class="sum-table"><caption>Top-10 genomes with the largest 16S rRNA copy numbers
</caption>
<tbody>
<tr><th class="alnleft">Organism</th><th class="numcol">Copy number</th><th class="alnleft">Assembly<br>accession</th></tr>
<tr>
<td colspan="4" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<div id="top-var" class="pad-anchor"></div>
<table class="sum-table"><caption>Top-10 genomes with the highest intragenomic variability of 16S rRNA genes
</caption>
<tbody>
<tr><th class="alnleft">Organism</th><th class="numcol">Sum of entropy<sup> *</sup> (bits)</th><th class="numcol">Mean entropy<sup> *</sup> (bits)</th><th class="numcol">Number of variable positions</th><th class="numcol">Gene copy number</th><th class="alnleft">Assembly<br>accession</th></tr>

<tr>
<td colspan="6" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_intragenomic_var_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['sum_entropy'] }}</td>
<td class="numcol">{{ row['mean_entropy'] }}</td>
<td class="numcol">{{ row['num_var_cols'] }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://www.ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="6" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_intragenomic_var_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['sum_entropy'] }}</td>
<td class="numcol">{{ row['mean_entropy'] }}</td>
<td class="numcol">{{ row['num_var_cols'] }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://www.ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<p><sup>*</sup> Entropy is <a href="https://en.wikipedia.org/wiki/Entropy_(information_theory)">Shannon entropy</a> calculated for each column of the multiple sequence alignment (MSA) of all full-length 16S rRNA genes of a genome. Entropy is then summed up (column “Sum of entropy”) and averaged (column “Mean entropy”).</p>
<p><sup>**</sup> <i>Halomicrobium</i> sp. ZPS1 is a quite remarkable case. This genome harbours two 16S rRNA genes, therefore entropy is equal to the number of mismatching nucleotides between sequences of the genes. Respectively, per cent of identity between these two gene sequences is 90.70%! This is remarkable because the usual (however arbitrary) <em>genus</em> demarcation threshold of per cent of identity is 95%.</p>

<div id="primers-coverages" class="pad-anchor"></div>
<table class="sum-table"><caption>Coverage<sup>*</sup> of primer pairs for different V-regions of bacterial 16S rRNA genes</caption>
<tbody class="primer-cov-tbody">
<tr>
  <th class="alnleft" rowspan="2">Phylum</th>
  <th class="numcol" rowspan="2">Number<br />of genomes</th>
  <th class="numcol">Full gene</th>
  <th class="numcol">V1–V2</th>
  <th class="numcol">V1–V3</th>
  <th class="numcol">V3–V4</th>
  <th class="numcol">V3–V5</th>
  <th class="numcol">V4</th>
  <th class="numcol">V4–V5</th>
  <th class="numcol">V4–V6</th>
  <th class="numcol">V5–V6</th>
  <th class="numcol">V5–V7</th>
  <th class="numcol">V6–V7</th>
  <th class="numcol">V6–V8</th>
</tr>
<tr>
  <th class="numcol">27F–1492R<br />(%)</th>
  <th class="numcol">27F–338R<br />(%)</th>
  <th class="numcol">27F–534R<br />(%)</th>
  <th class="numcol">341F–785R<br />(%)</th>
  <th class="numcol">341F–944R<br />(%)</th>
  <th class="numcol">515F–806R<br />(%)</th>
  <th class="numcol">515F–944R<br />(%)</th>
  <th class="numcol">515F–1100R<br />(%)</th>
  <th class="numcol">784F–1100R<br />(%)</th>
  <th class="numcol">784F–1193R<br />(%)</th>
  <th class="numcol">939F–1193R<br />(%)</th>
  <th class="numcol">939F–1378R<br />(%)</th>
</tr>
{% for _, row in ribogrove_primers_cov_df.iterrows() %}
<tr class="sumtab-row">
  <td><i>{{ row['Phylum'] }}</i></td>
  <td class="numcol">{{ row['num_genomes'] }}</td>
  <td class="numcol">{{ row['27F-1492R'] }}</td>
  <td class="numcol">{{ row['27F-338R'] }}</td>
  <td class="numcol">{{ row['27F-534R'] }}</td>
  <td class="numcol">{{ row['341F-785R'] }}</td>
  <td class="numcol">{{ row['341F-944R'] }}</td>
  <td class="numcol">{{ row['515F-806R'] }}</td>
  <td class="numcol">{{ row['515F-944R'] }}</td>
  <td class="numcol">{{ row['515F-1100R'] }}</td>
  <td class="numcol">{{ row['784F-1100R'] }}</td>
  <td class="numcol">{{ row['784F-1193R'] }}</td>
  <td class="numcol">{{ row['939F-1193R'] }}</td>
  <td class="numcol">{{ row['939F-1378R'] }}</td>
</tr>
{% endfor %}
</tbody>
</table>

<p><sup>*</sup> Coverage of a primer pair is the per cent of genomes having at least one 16S rRNA gene which can be amplified by PCR using this primer pair. For details, see our <a href="/cager/en/ribogrove#citing-ribogrove">paper about RiboGrove</a>.</p>
<p>You can find a more detailed table in the file <span class="samp">primer_pair_genomic_coverage.tsv</span> in the <a href="#downloads">metadata</a>. That table contains coverage not just for phyla, but also for each bacterial class, order, family, genus, and species. Moreover, that table contains coverage values for primer pair 1115F–1492R (V7–V9 region). In this table, it is omitted for brevity.</p>

<table class="sum-table"><caption>Primers used for coverage estimation</caption>
<tbody>
  <tr><th>Primer name</th><th>Sequence</th><th>Reference</th></tr>
  <tr class="sumtab-row"><td>27F</td><td>AGAGTTTGATYMTGGCTCAG</td><td><a href="https://doi.org/10.1128/AEM.02272-07">Frank et al., 2008</a></td></tr>
  <tr class="sumtab-row"><td>338R</td><td>GCTGCCTCCCGTAGGAGT</td><td><a href="https://doi.org/10.1128/aem.62.2.625-630.1996">Suzuki et al., 1996</a></td></tr>
  <tr class="sumtab-row"><td>341F<sup> *</sup></td><td>CCTACGGGNGGCWGCAG</td><td><a href="https://doi.org/10.1093/nar/gks808">Klindworth et al., 2013</a></td></tr>
  <tr class="sumtab-row"><td>515F</td><td>GTGCCAGCMGCCGCGGTAA</td><td><a href="https://doi.org/10.1111/j.1550-7408.1999.tb04612.x">Turner et al., 1999</a></td></tr>
  <tr class="sumtab-row"><td>534R</td><td>ATTACCGCGGCTGCTGG</td><td><a href="https://doi.org/10.1186/s40168-015-0087-4">Walker et al., 2015</a></td></tr>
  <tr class="sumtab-row"><td>784F</td><td>AGGATTAGATACCCTGGTA</td><td><a href="https://doi.org/10.1371/journal.pone.0002836">Andersson et al., 2008</a></td></tr>
  <tr class="sumtab-row"><td>785R<sup> *</sup></td><td>GACTACHVGGGTATCTAATCC</td><td><a href="https://doi.org/10.1093/nar/gks808">Klindworth et al., 2013</a></td></tr>
  <tr class="sumtab-row"><td>806R</td><td>GGACTACHVGGGTWTCTAAT</td><td><a href="https://doi.org/10.1073/pnas.1000080107">Caporaso et al., 2010</a></td></tr>
  <tr class="sumtab-row"><td>939F</td><td>GAATTGACGGGGGCCCGCACAAG</td><td><a href="https://doi.org/10.1016/j.anaerobe.2014.04.006">Lebuhn et al., 2014</a></td></tr>
  <tr class="sumtab-row"><td>944R</td><td>GAATTAAACCACATGCTC</td><td><a href="https://doi.org/10.1186/s40168-017-0396-x">Fuks et al., 2018</a></td></tr>
  <tr class="sumtab-row"><td>1100R</td><td>AGGGTTGCGCTCGTTG</td><td><a href="https://doi.org/10.1111/j.1550-7408.1999.tb04612.x">Turner et al., 1999</a></td></tr>
  <tr class="sumtab-row"><td>1193R</td><td>ACGTCATCCCCACCTTCC</td><td><a href="https://doi.org/10.1371/journal.pone.0056329">Bodenhausen et al, 2013</a></td></tr>
  <tr class="sumtab-row"><td>1378R</td><td>CGGTGTGTACAAGGCCCGGGAACG</td><td><a href="https://doi.org/10.1016/j.anaerobe.2014.04.006">Lebuhn et al., 2014</a></td></tr>
  <tr class="sumtab-row"><td>1492R</td><td>TACCTTGTTACGACTT</td><td><a href="https://doi.org/10.1128/AEM.02272-07">Frank et al., 2008</a></td></tr>
</tbody>
</table>

<p><sup>*</sup> Primers 341F and 785R are used in the <a href="https://support.illumina.com/downloads/16s_metagenomic_sequencing_library_preparation.html">protocol</a> for library preparation for sequencing of V3–V4 region of 16S rRNA genes on Illumina MiSeq.</p>

<hr>
<div id="searching-data" class="pad-anchor"></div>
<h2>Searching data in RiboGrove</h2>
<p>RiboGrove is a very minimalistic database — it comprises a collection of plain fasta files with metadata. Thus, extended search instruments are not available for it. We admit this problem and provide a list of suggestions below. The suggestions would help you to explore and select RiboGrove data.</p>
<div id="header-format" class="pad-anchor"></div>
<h3>Header format</h3>
<p>RiboGrove fasta data has the following format of header:</p>
<p class="samp-highl samp-vwide" style="font-size:14px">&gt;GCF_000978375.1:NZ_CP009686.1:8908-10459:plus ;d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__cereus; category:1</p>
<p>Major blocks of a header are separated by spaces. A header consists of three such blocks:</p>
<ol>
<li>Sequence ID (seqID): <span class="samp-highl">GCF_000978375.1:NZ_CP009686.1:8908-10459:plus</span>. SeqID, in turn, consists of four parts separated by semicolons (<span class="samp-highl">:</span>):
<ol>
<li>The Assembly accession of the genome from which the gene originates: <span class="samp-highl">GCF_000978375.1</span>.</li>
<li>The accession number of the RefSeq sequence, from which the gene originates: <span class="samp-highl">NZ_CP009686.1</span>.</li>
<li>Coordinates of the gene within this RefSeq genomic sequence: <span class="samp-highl">8908-10459</span> (coordinates are 1-based, left-closed and right-closed).</li>
<li>Strand of the RefSeq genomic sequence, where the gene is located: <span class="samp-highl">plus</span> (or <span class="samp-highl">minus</span>).</li>
</ol>
</li>
<li>A taxonomy string, comprising domain (<span class="samp-highl">Bacteria</span>), phylum (<span class="samp-highl">Firmicutes</span>), class (<span class="samp-highl">Bacilli</span>), order (<span class="samp-highl">Bacillales</span>), family (<span class="samp-highl">Bacillaceae</span>), genus (<span class="samp-highl">Bacillus</span>) names, and the specific epithet (<span class="samp-highl">cereus</span>).<br>
Each name is preceded by a prefix, which denotes rank: <span class="samp-highl">d__</span> for domain, <span class="samp-highl">p__</span> for phylum, <span class="samp-highl">c__</span> for class, <span class="samp-highl">o__</span> for order, <span class="samp-highl">f__</span> for family, <span class="samp-highl">g__</span> for genus, and <span class="samp-highl">s__</span> for specific epithet. Prefixes contain <strong>double</strong> underscores.<br>
The taxonomic names are separated and flanked by semicolons (<span class="samp-highl">;</span>).</li>
<li>The category of the genome, from which the gene sequence originates: (<span class="samp-highl">category:1</span>).</li>
</ol>
<h3>Sequence selection</h3>
<p>You can select specific sequences from fasta files using the <strong>Seqkit</strong> program (<a href="https://github.com/shenwei356/seqkit">GitHub repo</a>, <a href="https://bioinf.shenwei.me/seqkit/">documentation</a>). It is free, cross-platform, multifunctional and pretty fast and can process both gzipped and uncompressed fasta files. Programs <span class="samp">seqkit grep</span> and <span class="samp">seqkit seq</span> are useful for sequence selection.</p>
<div id="select-by-header" class="pad-anchor"></div>
<h4>Search sequences by header</h4>
<p>Given the downloaded fasta file <span class="samp">ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</span>, consider the following examples of sequence selection using <span class="samp">seqkit grep</span>:</p>
<p><strong>Example 1</strong>. Select a single sequence by SeqID.</p>
<p class="samp-highl samp-vwide">seqkit grep -p "GCF_000978375.1:NZ_CP009686.1:8908-10459:plus" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-p</span> option sets a pattern to search in fasta headers (only in sequence IDs, actually).</p>
<p><strong>Example 2</strong>. Select all gene sequences of a single RefSeq genomic sequence by accession number <span class="samp">NZ_CP009686.1</span>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ":NZ_CP009686.1:" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Here, two more options are required: <span class="samp">-n</span> and <span class="samp">-r</span>. The former tells the program to match the whole headers instead of IDs only. The latter tells the program to include partial matches into output, i.e. if the pattern is a substring of a header, the header will be printed to output.</p>
<p class="samp-comment">To ensure search specificity, surround the Accession.Version with colons (<span class="samp-highl">:</span>).</p>
<p><strong>Example 3</strong>. Select all gene sequences of a single genome (Assembly accession <a href="https://www.ncbi.nlm.nih.gov/assembly/GCF_019357495.1">GCF_019357495.1</a>).</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "GCF_019357495.1:" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">To ensure search specificity, put a colon (<span class="samp-highl">:</span>) after the assembly accession.</p>
<p><strong>Example 4</strong>. Select all actinobacterial sequences.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ";p__Actinobacteria;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">To ensure search specificity, surround the taxonomy name with semicolons (<span class="samp-highl">;</span>).</p>
<p><strong>Example 5</strong>. Select all sequences originating from category 1 genomes.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "category:1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Example 6</strong>. Select all sequences except for those belonging to <i>Firmicutes</i>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nvrp ";p__Firmicutes;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Recognize the <span class="samp">-v</span> option within the option sequence <span class="samp">-nvrp</span>. This option inverts match, i.e. output will comprise sequences, headers of which do <strong>not</strong> contain the substring “<span class="samp">;p__Firmicutes;</span>”.</p>
<div id="select-by-len" class="pad-anchor"></div>
<h4>Search sequences by length</h4>
<p>You can use the <span class="samp">seqkit seq</span> program to select sequences by length.</p>
<p><strong>Example 1</strong>. Select all sequences longer than 1600 bp.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1601 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-m</span> option sets the minimum length of a sequence to be printed to output.</p>
<p><strong>Example 2</strong>. Select all sequences shorter than 1500 bp.</p>
<p class="samp-highl samp-vwide">seqkit seq -M 1499 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-M</span> option sets the maximum length of a sequence to be printed to output.</p>
<p><strong>Example 3</strong>. Select all sequences having length in range [1500, 1600] bp.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1500 -M 1600 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<div id="select-head" class="pad-anchor"></div>
<h4>Selecting header data</h4>
<p>It is sometimes useful to retrieve only header information from a fasta file. You can use the <span class="samp">seqkit seq</span> program for it.</p>
<p><strong>Example 1</strong>. Select all headers.</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-n</span> option tells the program to output only headers.</p>
<p><strong>Example 2</strong>. Select all SeqIDs (header parts before the first space).</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-i</span> option tells the program to output only sequence IDs.</p>
<p><strong>Example 3</strong>. Select all RefSeq “Assession.Version”s.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f2 -d':' | sort | uniq</p>
<p class="samp-comment">This might be done only if you have <span class="samp">cut</span>, <span class="samp">sort</span> and <span class="samp">uniq</span> utilities installed (Linux and Mac OS systems should have them built-in).</p>
<p><strong>Example 4</strong>. Select all Assembly accessions.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f1 -d':' | sort | uniq</p>
<p class="samp-comment">This might be done only if you have <span class="samp">cut</span>, <span class="samp">sort</span> and <span class="samp">uniq</span> utilities installed (Linux and Mac OS systems should have them built-in).</p>
<p><strong>Example 5</strong>. Select all phylum names.</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | grep -Eo ';p__[^;]+' | sed -E 's/;|p__//g' | sort | uniq</p>
<p class="samp-comment">This might be done only if you have <span class="samp">grep</span>, <span class="samp">sed</span>, <span class="samp">sort</span> and <span class="samp">uniq</span> utilities installed (Linux and Mac OS systems should have them built-in).</p>
<hr>
<div id="contacts" class="pad-anchor"></div>
<h2>Contacts</h2>
<p>For any questions concerning RiboGrove, please contact Maxim Sikolenko at sikolenko<img class="sabaka" src="/wp-content/uploads/cager/pes2.gif" alt="[ at ]" align="bottom">bio.bsu.by.</p>
<hr>
<div id="citing-ribogrove" class="pad-anchor"></div>
<h2>Citing RiboGrove</h2>
<p>If you find RiboGrove useful for your research please cite:</p>
<p>Maxim A. Sikolenko, Leonid N. Valentovich. “RiboGrove: a database of full-length prokaryotic 16S rRNA genes derived from completely assembled genomes” // Research in Microbiology, Volume 173, Issue 4, May 2022, 103936.<br> (DOI: <a href="https://doi.org/10.1016/j.resmic.2022.103936">10.1016/j.resmic.2022.103936</a>).</p>
<hr>
<div id="faq" class="pad-anchor"></div>
<h2>Questions people ask about RiboGrove</h2>
<div><details>
  <summary><b>1. How do I create QIIME2-compatible taxonomy file from RiboGrove data?</b></summary>
  <p>Please use the <a href="https://github.com/masikol/ribogrove-tools/blob/main/data_transformation_scripts/make_qiime_taxonomy_file.py"><span class="samp">make_qiime_taxonomy_file.py</span> script</a> to convert the RiboGrove file <span class="samp">metadata/taxonomy.tsv</span> to a QIIME2-compatible file. You can find out how to use this script in the corresponding <a href="https://github.com/masikol/ribogrove-tools/blob/main/data_transformation_scripts/README.md">README file</a>.</p>
</details></div><br>
<div><details>
  <summary><b>2. How do I save selected sequences in Seqkit to a file?</b></summary>
  <p>People have already provided several useful answers in the corresponding discussion: <a href="https://bioinformatics.stackexchange.com/questions/20915/how-do-i-save-selected-sequences-in-seqkit-to-a-file">https://bioinformatics.stackexchange.com/questions/20915/how-do-i-save-selected-sequences-in-seqkit-to-a-file</a>.</p>
</details></div><br>
<div><details>
  <summary><b>3. How do I search a FASTA database by sequence in Seqkit?</b></summary>
  <p>People have already provided several useful answers in the corresponding discussion: <a href="https://www.biostars.org/p/9561418">https://www.biostars.org/p/9561418</a>.</p>
</details></div>
