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
<li><a href="#citing-ribogrove">Citing RiboGrove</a></li>
<li><a href="#downloads">Downloads</a>
<ul class="ribogrove-nested-list">
  <li><a href="#latest-release">Latest RiboGrove release</li>
  <li><a href="#release-archive">RiboGrove release archive</li>
</ul>
</li>
<li><a href="#summary">Statistical summary</a>
<ul class="ribogrove-nested-list">
<li><a href="#db-size">RiboGrove size</a></li>
<li><a href="#gene-lens">16S rRNA gene lengths</a></li>
<li><a href="#copy-number">16S rRNA gene copy number</a></li>
<li><a href="#top-longest">Top-10 longest 16S rRNA genes</a></li>
<li><a href="#top-shortest">Top-10 shortest 16S rRNA genes</a></li>
<li><a href="#top-copy-num">Top-10 genomes with largest 16S rRNA copy numbers</a></li>
<li><a href="#top-var">Top-10 genomes with highest intragenomic variability of 16S rRNA genes</a></li>
</ul>
</li>
<li><a href="#searching-data">Searching data in RiboGrove</a>
<ul class="ribogrove-nested-list">
<li><a href="#header-format">Header format</a></li>
<li><a href="#subset-by-header">Subsetting sequences by header</a></li>
<li><a href="#subset-by-len">Subsetting sequences by length</a></li>
<li><a href="#select-head">Selecting header data</a></li>
</ul>
</li>
</ul>
<hr>
<div id="overview" class="pad-anchor"></div>
<h2>What is RiboGrove</h2>
<p>RiboGrove is a database of 16S rRNA gene sequences of bacteria and archaea.</p>
<p>RiboGrove is based on the <a href="https://ncbi.nlm.nih.gov/refseq/">RefSeq</a> database. It contains only full-length sequences of 16S rRNA genes, and the sequences are derived from completely assembled prokaryotic genomes deposited in RefSeq. Hence we posit high reliability of RiboGrove sequences.</p>
<div id="other-databases" class="pad-anchor"></div>
<h3>RiboGrove and other 16S rRNA databases</h3>
<p>Here is a summary showing what is the (qualitative) difference between RiboGrove and similar rRNA sequence databases, namely <a href="https://www.arb-silva.de/">Silva</a>, <a href="http://rdp.cme.msu.edu/">RDP</a>, and <a href="https://greengenes.secondgenome.com/">Greengenes</a>. Briefly, RiboGrove is inferior in sequence amount and diversity, but superior in sequence reliability.</p>
<table class="sum-table">
<tbody>
  <tr class="sumtab-row"><th></th><th class="alcentr">RiboGrove</th><th class="alcentr">Silva</th><th class="alcentr">RDP</th><th class="alcentr">Greengenes</th></tr>
  <tr class="sumtab-row"><td class="feature-col">Represented organisms</td><td class="alcentr">Bacteria<br>Archaea</td><td class="alcentr">Bacteria<br>Archaea<br>Eukaryotes</td><td class="alcentr">Bacteria<br>Archaea<br>Eukaryotes</td><td class="alcentr">Bacteria<br>Archaea</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Represented ribosome subunits</td><td class="alcentr">Small</td><td class="alcentr">Large<br>Small</td><td class="alcentr">Large<br>Small</td><td class="alcentr">Small</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Contains sequences from assembles genomes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Contains amplicon sequences</td><td class="alcentr">No</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Contains partial gene sequences</td><td class="alcentr">No</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td><td class="alcentr">Yes</td></tr>
</tbody>
</table>
<div id="genome-categories" class="pad-anchor"></div>
<h3>Genome categories</h3>
<p>All genomes used for RiboGrove construction were divided into three categories according to their expected reliability:</p>
<ol style="list-style:none;">
  <li><strong>Category 1</strong> (the highest reliability). Genomes showing no signs of a low-quality assembly and sequenced either with PacBio technology or with combination (Oxford&nbsp;Nanopore&nbsp;+&nbsp;Illumina).</li>
  <li><strong>Category 2.</strong> Genomes showing no signs of a low-quality assembly and sequenced with any other technology (or the technology is not specified).</li>
  <li><strong>Category 3</strong> (the lowest reliability). Genomes showing at least one sign of a low-quality assembly.</li>
</ol>
<p>Signs of a low-quality assembly are the following:</p>
<ul>
  <li>The genome contains degenerate base(s) in 16S rRNA genes sequences.</li>
  <li>The assembly includes at least one RefSeq record whose title contains the phrase <q>map unlocalized</q> and this record contains a 16S rRNA gene or a part of it.</li>
</ul>
<br>
<p>The software used for the RiboGrove construction can be found in the following GitHub repository: <a href="https://github.com/masikol/ribogrove-tools">ribogrove-tools</a>.</p>
<hr>
<div id="citing-ribogrove" class="pad-anchor"></div>
<h2>Citing RiboGrove</h2>
<p>If you find RiboGrove useful for your research please cite:</p>
<p>Maxim A. Sikolenko, Leonid N. Valentovich. &ldquo;RiboGrove: a database of full-length prokaryotic 16S rRNA genes derived from completely assembled genomes&rdquo; // Research in Microbiology, Volume 173, Issue 4, May 2022, 103936.<br> (DOI: <a href="https://doi.org/10.1016/j.resmic.2022.103936">10.1016/j.resmic.2022.103936</a>).</p>
<hr>
<div id="downloads" class="pad-anchor"></div>
<h2>Downloads</h2>
<div id="latest-release" class="pad-anchor"></div>
<h3>Latest RiboGrove release &mdash; {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
<p>The release is based on RefSeq release {{ refseq_release }}.</p>
<ul>
<li>A fasta file of full-length 16S gene sequences. <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Download (gzipped fasta file, {{ final_fasta_fsize_fmt }} MB)</a></li>
<li>A &ldquo;raw&ldquo; version of the fasta file above. This file contains partial sequences. <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/raw_ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Download (gzipped fasta file, {{ raw_fasta_fsize_fmt }} MB)</a></li>
<li>Metadata <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/metadata_ribogrove_{{ ribogrove_release_number }}.zip">Download (zip archive {{ metadata_fsize_fmt }} MB)</a><br>
  <details>
    <summary><small>What information exactly does the metadata contain?</small></summary>
    <small>
      <p>The metadata consists of the following files:</p>
      <ol>
      <li>
        <span class="samp-highl">source_RefSeq_genomes.tsv</span><br>
        This is a TSV files, which contains information about what genomes were used for the RiboGrove construction.
      </li>
      <li>
        <span class="samp-highl">gene_seqs_statistics.tsv</span>, <span class="samp-highl">raw_gene_seqs_statistics.tsv</span><br>
        These are TSV file, which contains nucleotide conposition, size, genomic and taxonomic affiliation of the gene sequences. The &ldquo;raw&rdquo; file also includes information about partial genes.
      </li>
      <li>
        <span class="samp-highl">per_replicon_statistics.tsv</span>, <span class="samp-highl">raw_per_replicon_statistics.tsv</span><br>
        These are TSV files, which contain information about number of 16S rRNA genes in each RefSeq genomic sequences, and about sizes of these genes. The &ldquo;raw&rdquo; file also includes information about partial genes.
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
        <span class="samp-highl">cmscan_output_table.tblout</span><br>
        This is a TSV file, which contains complete output of <a href="http://eddylab.org/infernal/">cmscan</a> program outputted after processing all 16S rRNA sequences.
      </li>
      <li>
        <span class="samp-highl">entropy_summary.tsv</span><br>
        This is a TSV file, which contains summary of instragenomic variability of the 16S rRNA genes. Intragenomic variability are calculated only for the category&nbsp;1 genomes having more than one 16S rRNA gene. Intragenomic variability is evaluated using Shannon entropy. We align gene sequences from each genome using <a href="https://drive5.com/muscle5/">MUSCLE</a>, and then we calculate Shannon entropy for each multiple alignment column (i.e. base).
      </li>
      <li>
        <span class="samp-highl">QIIME2-compatible-taxonomy.txt</span><br>
        This is a TSV file, which can be used to train a QIIME2 classifier (see the <a href="https://docs.qiime2.org">tutorial</a>).
      </li>
    </ol></small>
  </details>
</li>
</ul>
<p>The fasta files are compressed with gzip, and the metadata file is a zip archive. To uncompress them, Linux and Mac OS users may use gzip and zip programs, they should be built-in. For Windows users, the free and open source (de)compression program <a href="https://www.7-zip.org/">7-Zip</a> is available.</p>
<div id="release-archive" class="pad-anchor"></div>
<h3>RiboGrove release archive</h3>
<p>You can find all releases at the <a href="/cager/en/content/60-ribogrove-archive">RiboGrove release archive</a>.</p>
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
<p><sup>*</sup> Metrics marked with this sign were calculated with preliminary normalization, i.e. median within-species gene length were used for the summary.</p>
<div id="copy-number" class="pad-anchor"></div>
<table class="sum-table">
<caption>16S rRNA gene copy number<br><small>(Bacteria and Archaea)</small></caption>
<tbody>
<tr><th class="numcol">Copy number <sup>*</sup></th><th class="numcol">Number of species</th><th class="numcol">Per cent of species (%)</th></tr>

{% for _, row in ribogrove_copy_number_df.iterrows() %}
<tr class="sumtab-row">
<td class="numcol">{{ row['copy_number'] }}</td>
<td class="numcol">{{ row['number_of_species'] }}</td>
<td class="numcol">{{ row['percent_of_species'] }}</td>
</tr>
{% endfor %}

</tbody>
</table>
<p><sup>*</sup> These are median within-species copy numbers.</p>

<div id="top-longest" class="pad-anchor"></div>
<table class="sum-table"><caption>Top-10 longest 16S rRNA genes
</caption>
<tbody>
<tr><th class="alnleft">Organism</th><th class="numcol">Gene length (bp)</th><th class="alnleft">RiboGrove Sequence ID(s)</th><th class="alnleft">Assembly ID</th></tr>
<tr>
<td colspan="4" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_longest_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['title']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_longest_df.query('domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['title']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<div id="top-shortest" class="pad-anchor"></div>
<table class="sum-table"><caption>Top-10 shortest 16S rRNA genes
</caption>
<tbody>
<tr><th class="alnleft">Organism</th><th class="numcol">Gene length (bp)</th><th class="alnleft">RiboGrove Sequence ID(s)</th><th class="alnleft">Assembly ID</th></tr>
<tr>
<td colspan="4" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_shortest_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['title']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_shortest_df.query('domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['title']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<div id="top-copy-num" class="pad-anchor"></div>
<table class="sum-table"><caption>Top-10 genomes with largest 16S rRNA copy numbers
</caption>
<tbody>
<tr><th class="alnleft">Organism</th><th class="numcol">Copy number</th><th class="alnleft">Assembly ID</th></tr>
<tr>
<td colspan="4" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['title']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['title']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<div id="top-var" class="pad-anchor"></div>
<table class="sum-table"><caption>Top-10 genomes with highest intragenomic variability of 16S rRNA genes
</caption>
<tbody>
<tr><th class="alnleft">Organism</th><th class="numcol">Sum of entropy<sup> *</sup> (bits)</th><th class="numcol">Mean entropy<sup> *</sup> (bits)</th><th class="numcol">Number of variable positions</th><th class="numcol">Gene copy number</th><th class="alnleft">Assembly ID</th></tr>

<tr>
<td colspan="6" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_intragenomic_var_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['title']) }}</td>
<td class="numcol">{{ row['sum_entropy'] }}</td>
<td class="numcol">{{ row['mean_entropy'] }}</td>
<td class="numcol">{{ row['num_var_cols'] }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://www.ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="6" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_intragenomic_var_df.query('domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['title']) }}</td>
<td class="numcol">{{ row['sum_entropy'] }}</td>
<td class="numcol">{{ row['mean_entropy'] }}</td>
<td class="numcol">{{ row['num_var_cols'] }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://www.ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<p><sup>*</sup> Entropy is <a href="https://en.wikipedia.org/wiki/Entropy_(information_theory)">Shannon entropy</a> calculated for each column of the multiple sequence alignment (MSA) of all full-length 16S rRNA genes of a genome. Entropy is then summed up (column &ldquo;Sum of entropy&rdquo;) and averaged (column &ldquo;Mean entropy&rdquo;).</p>
<p><sup>**</sup> <i>Halomicrobium</i> sp. ZPS1 is a quite remarkable case. This genome harbours two 16S rRNA genes, therefore entropy is equal to the number of mismatching nucleotides between sequences of the genes. Respectively, per cent of identity between these two gene sequences is 90.70%! This is remarkable, because the usual (however arbitrary) <em>genus</em> demarcation threshold of per cent of identity is 95%.</p>

<hr>
<div id="searching-data" class="pad-anchor"></div>
<h2>Searching data in RiboGrove</h2>
<p>RiboGrove is a very minimalistic database &mdash; it comprises a collection of plain fasta files with metadata. Thus, extended search instruments are not available for it. We admit this problem and provide a list of suggestions below. The suggestions would help you to explore and subset RiboGrove data.</p>
<div id="header-format" class="pad-anchor"></div>
<h3>Header format</h3>
<p>RiboGrove fasta data has the following format of header:</p>
<p class="samp-highl samp-vwide">&gt;NZ_CP079719.1:86193-87742_plus Bacillus_velezensis ;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus; category:2</p>
<p>Major blocks of a header are separated by spaces. A header consists of four such blocks:</p>
<ol>
<li>Sequence ID (seqID): <span class="samp-highl">NZ_CP079719.1:86193-87742_plus</span>. SeqID, in turn, consists of three parts:<ol>
<li>the accession number of the RefSeq sequence, from which the gene originates: <span class="samp-highl">NZ_CP079719.1</span>;</li>
<li>coordinates of the gene within this RefSeq genomic sequence: <span class="samp-highl">86193-87742</span> (coordinates are 1-based, left-closed and right-closed);</li>
<li>strand of the RefSeq genomic sequence, where the gene is located: <span class="samp-highl">plus</span> (or <span class="samp-highl">minus</span>).</li>
</ol></li>
<li>The complete name of the organism, according to the <a href="https://www.ncbi.nlm.nih.gov/taxonomy">NCBI Taxonomy database</a>: <span class="samp-highl">Bacillus_velezensis</span>.</li>
<li>Taxonomy string, comprising domain (<span class="samp-highl">Bacteria</span>), phylum (<span class="samp-highl">Firmicutes</span>), class (<span class="samp-highl">Bacilli</span>), order (<span class="samp-highl">Bacillales</span>), family (<span class="samp-highl">Bacillaceae</span>), and genus (<span class="samp-highl">Bacillus</span>) names. The names are separated and flanked by semicolons (<span class="samp-highl">;</span>).</li>
<li>The category of the genome, from which the gene sequence originates: (<span class="samp-highl">category:2</span>. Category 1 genomes are of the highest reliability, and category 3 genomes are least reliable.</li>
</ol>
<h3>Subsetting sequences</h3>
<p>You can subset specific sequences from fasta files using the <strong>Seqkit</strong> program (<a href="https://github.com/shenwei356/seqkit">GitHub repo</a>, <a href="https://bioinf.shenwei.me/seqkit/">documentation</a>). It is free, cross-platform, multifunctional and pretty fast and can process both gzipped and uncompressed fasta files. Programs <span class="samp">seqkit grep</span> and <span class="samp">seqkit seq</span> are useful for sequence subsetting.</p>
<div id="subset-by-header" class="pad-anchor"></div>
<h4>Subsetting sequences by header</h4>
<p>Given the downloaded fasta file <span class="samp">ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</span>, consider the following examples of sequence subsetting using <span class="samp">seqkit grep</span>:</p>
<p><strong>Example 1</strong>. Select a single sequence by SeqID.</p>
<p class="samp-highl samp-vwide">seqkit grep -p "NZ_CP079719.1:86193-87742_plus" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-p</span> option sets a pattern to search in fasta headers (only in sequence IDs, actually).</p>
<p><strong>Example 2</strong>. Select all gene sequences of a single RefSeq genomic sequence by accession number.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "NZ_CP079719.1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Here, two more options are required: <span class="samp">-n</span> and <span class="samp">-r</span>. The former tells the program to match the whole headers instead of IDs only. The latter tells the program not to exclude partial matches from output, i.e. if the pattern is a substring of a header, the header will be printed to output.</p>
<p><strong>Example 3</strong>. Select all gene sequences of a single genome (Assembly ID <a href="https://www.ncbi.nlm.nih.gov/assembly/10577151">10577151</a>), which has two replicons: <a href="https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP079110.1">NZ_CP079110.1</a> and <a href="https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP079111.1">NZ_CP079111.1</a>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nr -p "NZ_CP079110.1" -p "NZ_CP079111.1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Example 4</strong>. Select all actinobacterial sequences.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ";Actinobacteria;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Just in case, surround the taxonomy name with semicolons (<span class="samp-highl">;</span>).</p>
<p><strong>Example 5</strong>. Select all sequences originating from category 1 genomes.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "category:1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Example 6</strong>. Select all sequences except for those belonging to <i>Firmicutes</i>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nvrp ";Firmicutes;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Recognize the <span class="samp">-v</span> option within option sequence <span class="samp">-nvrp</span>. This option inverts match, i.e. without it the search would result in sequences belonging to <i>Firmicutes</i> only.</p>
<div id="subset-by-len" class="pad-anchor"></div>
<h4>Subsetting sequences by length</h4>
<p>You can use program <span class="samp">seqkit seq</span> to subset sequences by length.</p>
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
<p>It is sometimes useful to retrieve only header information from a fasta file. You can use program <span class="samp">seqkit seq</span> for it.</p>
<p><strong>Example 1</strong>. Select all headers.</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-n</span> option tells the program to output only headers.</p>
<p><strong>Example 2</strong>. Select all SeqIDs (header parts before the first space).</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-i</span> option tells the program to output only sequence IDs.</p>
<p><strong>Example 3</strong>. Select all accession numbers.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f1 -d':' | sort | uniq</p>
<p class="samp-comment">This might be done only if you have <span class="samp">cut</span>, <span class="samp">sort</span> and <span class="samp">uniq</span> utilities installed (Linux and Mac OS systems should have them built-in).</p>
