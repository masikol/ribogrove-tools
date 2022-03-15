<img src="/cager/img/RiboGrove_logo.png"><br /><br />
<h2>Contents</h2>
<ul>
<li><a href="#overview">Overview</a></li>
<li><a href="#downloads">Downloads</a></li>
<li><a href="#summary">Summary</a>
<ul style="list-style-type: none;">
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
<ul style="list-style-type: none;">
<li><a href="#subset-by-header">Subsetting sequences by header</a></li>
<li><a href="#subset-by-len">Subsetting sequences by length</a></li>
<li><a href="#select-head">Selecting header data</a></li>
</ul>
</li>
<li><a href="#citing-ribogrove">Citing RiboGrove</a></li>
</ul>
<hr />
<h2>Overview</h2>
<div id="overview" class="pad-anchor"></div>
<p>RiboGrove is a database of 16S rRNA gene sequences of bacteria and archaea. RiboGrove contains only full-length sequences of 16S rRNA genes since amplicon data is not used for its construction. In other words, unlike similar rRNA resources like <a href="https://www.arb-silva.de/">Silva</a> and <a href="http://rdp.cme.msu.edu/">RDP</a>, RiboGrove does not contain partial gene sequences.</p>
<p>Completely assembled prokaryotic genomes deposited in <a href="https://www.ncbi.nlm.nih.gov/refseq/">RefSeq</a> database were used for the RiboGrove construction.</p>
<p>The software used for the RiboGrove construction can be found in the following GitHub repository: <a href="https://github.com/masikol/ribogrove-paper">ribogrove-paper</a>.</p>
<hr />
<h2>Downloads</h2>
<div id="downloads" class="pad-anchor"></div>
<h3>Latest release &mdash; {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
<p>The release is based on RefSeq release {{ refseq_release }}.</p>
<ul>
<li>A fasta file of full-length 16S gene sequences. <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Download (gzipped fasta file, {{ final_fasta_fsize_fmt }} MB)</a></li>
<li>A "raw" version of the fasta file above. This file contains partial sequences. <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/raw_ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Download (gzipped fasta file, {{ raw_fasta_fsize_fmt }} MB)</a></li>
<li>Metadata <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/metadata_ribogrove_{{ ribogrove_release_number }}.zip">Download (zip archive {{ metadata_fsize_fmt }} MB)</a></li>
</ul>
<p>The fasta files are compressed with gzip, and the metadata file if a zip archive. To uncompress them, Linux and Mac OS users may use gzip and zip programs, they should be built-in. For Windows users, the free and open source (de)compression program <a href="https://www.7-zip.org/">7-Zip</a> is available.</p>
<h3>RiboGrove release archive</h3>
<p>You can find all releases at the <a href="/cager/en/content/60-ribogrove-archive">RiboGrove release archive</a>.</p>
<hr />
<h2>Summary</h2>
<div id="summary" class="pad-anchor"></div>
<table class="sum-table"><caption>RiboGrove size
<div id="db-size" class="pad-anchor"></div>
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
<table class="sum-table"><caption>16S rRNA gene lengths
<div id="gene-lens" class="pad-anchor"></div>
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
<p><sup>*</sup> &mdash; metrics marked with this sign were calculated with preliminary normalization, i.e. median within-species gene length were used for the summary.</p>
<table class="sum-table"><caption>16S rRNA gene copy number
<div id="copy-number" class="pad-anchor"></div>
</caption>
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
<p><sup>*</sup> — these are median within-species copy numbers.</p>

<table class="sum-table"><caption>Top-10 longest 16S rRNA genes
<div id="top-longest" class="pad-anchor"></div>
</caption>
<tbody>
<tr><th class="alnleft">Assembly ID</th><th class="numcol">Gene length (bp)</th><th class="alnleft">RiboGrove Sequence ID(s)</th><th class="alnleft">Organism</th></tr>
<tr>
<td colspan="4" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_longest_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br />'.join(row['seqID']) }}</td>
<td><i>{{ row['species'] }}</i></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_longest_df.query('domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br />'.join(row['seqID']) }}</td>
<td><i>{{ row['species'] }}</i></td>
</tr>
{% endfor %}
</tbody>
</table>

<table class="sum-table"><caption>Top-10 shortest 16S rRNA genes
<div id="top-shortest" class="pad-anchor"></div>
</caption>
<tbody>
<tr><th class="alnleft">Assembly ID</th><th class="numcol">Gene length (bp)</th><th class="alnleft">RiboGrove Sequence ID(s)</th><th class="alnleft">Organism</th></tr>
<tr>
<td colspan="4" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_shortest_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br />'.join(row['seqID']) }}</td>
<td><i>{{ row['species'] }}</i></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_shortest_df.query('domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br />'.join(row['seqID']) }}</td>
<td><i>{{ row['species'] }}</i></td>
</tr>
{% endfor %}
</tbody>
</table>

<table class="sum-table"><caption>Top-10 genomes with largest 16S rRNA copy numbers
<div id="top-copy-num" class="pad-anchor"></div>
</caption>
<tbody>
<tr><th class="alnleft">Assembly ID</th><th class="numcol">Copy number</th><th class="alnleft">Organism</th></tr>

<tr>
<td colspan="4" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><i>{{ row['species'] }}</i></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><i>{{ row['species'] }}</i></td>
</tr>
{% endfor %}
</tbody>
</table>

<table class="sum-table"><caption>Top-10 genomes with highest intragenomic variability of 16S rRNA genes
<div id="top-var" class="pad-anchor"></div>
</caption>
<tbody>
<tr><th class="alnleft">Assembly ID</th><th class="numcol">Sum of entropy<sup>*</sup> (bits)</th><th class="numcol">Mean entropy<sup>*</sup> (bits)</th><th class="numcol">Number of variable positions</th><th class="numcol">Copy number</th><th class="alnleft">Organism</th></tr>

<tr>
<td colspan="6" class="subhead">Bacteria</td>
</tr>
{% for _, row in ribogrove_top_intragenomic_var_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td><a href="https://www.ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
<td class="numcol">{{ row['sum_entropy'] }}</td>
<td class="numcol">{{ row['mean_entropy'] }}</td>
<td class="numcol">{{ row['num_var_cols'] }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><i>{{ row['species'] }}</i></td>
</tr>
{% endfor %}
<tr>
<td colspan="6" class="subhead">Archaea</td>
</tr>
{% for _, row in ribogrove_top_intragenomic_var_df.query('domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td><a href="https://www.ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
<td class="numcol">{{ row['sum_entropy'] }}</td>
<td class="numcol">{{ row['mean_entropy'] }}</td>
<td class="numcol">{{ row['num_var_cols'] }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><i>{{ row['species'] }}</i></td>
</tr>
{% endfor %}
</tbody>
</table>

<p><sup>*</sup> — entropy is <a href="https://en.wikipedia.org/wiki/Entropy_(information_theory)">Shannon entropy</a> calculated for each column of the multiple sequence alignment (MSA) of all full-length 16S rRNA genes of a genome. Entropy is then summed up (column "Sum of entropy") and averaged (column "Mean entropy").</p>
<hr />
<h2>Searching data in RiboGrove</h2>
<div id="searching-data" class="pad-anchor"></div>
<p>RiboGrove is a very minimalistic database — it comprises a collection of plain fasta files with metadata. Thus, extended search instruments are not available for it. We admit this problem and provide a list of suggestions below. The suggestions would help you to explore and subset RiboGrove data.</p>
<h3>Header format</h3>
<p>RiboGrove fasta data has the following format of header:</p>
<p class="samp-highl samp-vwide">&gt;NZ_CP079719.1:86193-87742_plus Bacillus_velezensis ;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus; category:2</p>
<p>Major blocks of a header are separated by spaces. A header consists of four such blocks:</p>
<ol>
<li>Sequence ID (seqID): <span class="samp-highl">NZ_CP079719.1:86193-87742_plus</span>. SeqID, in turn, consists of three parts:<ol>
<li>the accession number of the RefSeq sequence, from which the gene originates: <span class="samp-highl">NZ_CP079719.1</span>;</li>
<li>coordinates of the gene within this RefSeq genomic sequence: <span class="samp-highl">86193-87742</span> (coordinates are 1-based, left-closed and right-closed);</li>
<li>strand of the RefSeq genomic sequence, where the gene is located: <span class="samp-highl">plus</span>.</li>
</ol></li>
<li>The complete name of the organism, according to the <a href="https://www.ncbi.nlm.nih.gov/taxonomy">NCBI Taxonomy database</a>: <span class="samp-highl">Bacillus_velezensis</span>.</li>
<li>Taxonomy string, comprising domain (<span class="samp-highl">Bacteria</span>), phylum (<span class="samp-highl">Firmicutes</span>), class (<span class="samp-highl">Bacilli</span>), order (<span class="samp-highl">Bacillales</span>), family (<span class="samp-highl">Bacillaceae</span>), and genus (<span class="samp-highl">Bacillus</span>) names. The names are separated and flanked by semicolons (<span class="samp-highl">;</span>).</li>
<li>The category of the genome, from which the gene sequence originates: (<span class="samp-highl">category:2</span>. Category 1 genomes are of the highest reliability, and category 3 genomes are least reliable.</li>
</ol>
<h3>Subsetting sequences</h3>
<p>You can subset specific sequences from fasta files using the <strong>Seqkit</strong> program (<a href="https://github.com/shenwei356/seqkit">GitHub repo</a>, <a href="https://bioinf.shenwei.me/seqkit/">documentation</a>). It is free, cross-platform, multifunctional and pretty fast and can process both gzipped and uncompressed fasta files. To be specific, the <span class="samp">seqkit grep</span> is the program that should be used for subsetting.</p>
<h4>Subsetting sequences by header</h4>
<div id="subset-by-header" class="pad-anchor"></div>
<p>Given the downloaded fasta file <span class="samp">ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</span>, consider the following examples of sequence subsetting using <span class="samp">seqkit grep</span>:</p>
<p><strong>Example 1</strong>. Select a single sequence by SeqID</p>
<p class="samp-highl samp-vwide">seqkit grep -p "NZ_CP079719.1:86193-87742_plus" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span>-p</span> option sets a pattern to search in fasta headers (only in sequence IDs, actually).</p>
<p><strong>Example 2</strong>. Select all gene sequences of a single RefSeq genomic sequence by accession number</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "NZ_CP079719.1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Here, two more options are required: <span class="samp">-n</span> and <span class="samp">-r</span>. The former tells the program to match the whole headers instead of IDs only. The latter tells the program not to exclude partial matches from output, i.e. if the pattern is a substring of a header, the header will be printed to output.</p>
<p><strong>Example 3</strong>. Select all gene sequences of a single genome (Assembly ID <a href="https://www.ncbi.nlm.nih.gov/assembly/10577151">10577151</a>), which has two replicons: <a href="https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP079110.1">NZ_CP079110.1</a> and <a href="https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP079111.1">NZ_CP079111.1</a></p>
<p class="samp-highl samp-vwide">seqkit grep -nr -p "NZ_CP079110.1" -p "NZ_CP079111.1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Example 4</strong>. Select all actinobacterial sequences</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ";Actinobacteria;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Just in case, surround the taxonomy name with semicolons (<span class="samp-highl">;</span>).</p>
<p><strong>Example 5</strong>. Select all sequences originating from category 1 genomes</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "category:1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Example 6</strong>. Select all sequences except for those belonging to <i>Firmicutes</i></p>
<p class="samp-highl samp-vwide">seqkit grep -nvrp ";Firmicutes;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Recognize the <span class="samp">-v</span> option within option sequence <span class="samp">-nvrp</span>. This option inverts match, i.e. without it the search would result in sequences belonging to <i>Firmicutes</i> only.</p>
<h4>Subsetting sequences by length</h4>
<div id="subset-by-len" class="pad-anchor"></div>
<p>Given the downloaded fasta file <span class="samp">ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</span>, consider the following examples of sequence subsetting using <span class="samp">seqkit seq</span>:</p>
<p><strong>Example 1</strong>. Select all sequences longer than 1600 bp</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1601 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-m</span> option sets the minimum length of a sequence to be printed to output.</p>
<p><strong>Example 2</strong>. Select all sequences shorter than 1500 bp</p>
<p class="samp-highl samp-vwide">seqkit seq -M 1499 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-M</span> option sets the maximum length of a sequence to be printed to output.</p>
<p><strong>Example 3</strong>. Select all sequences having length in range [1500, 1600] bp</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1500 -M 1600 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<h4>Selecting header data</h4>
<div id="select-head" class="pad-anchor"></div>
<p>It is sometimes useful to retrieve only header information from a fasta file. So, please consider the following examples of data selection using <span class="samp">seqkit seq</span>:</p>
<p><strong>Example 1</strong>. Select all headers</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-n</span> option tells the program to output only headers.</p>
<p><strong>Example 2</strong>. Select all SeqIDs (header parts before the first space)</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">The <span class="samp">-i</span> option tells the program to output only sequence IDs.</p>
<p><strong>Example 3</strong>. Select all accession numbers</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | sort | uniq</p>
<p class="samp-comment">This might be done only if you have <span class="samp">sort</span> and <span class="samp">uniq</span> utilities installed (Linux and Mac OS systems should have them built-in).</p>
<h2>Citing RiboGrove</h2>
<div id="citing-ribogrove" class="pad-anchor"></div>
<p>If you find RiboGrove useful for your research please cite:</p>
<p>Maxim A. Sikolenko, Leonid N. Valentovich. "RiboGrove: a database of full-length prokaryotic 16S rRNA genes derived from completely assembled genomes" // Research in Microbiology, 2022, 103936.<br /> (DOI: <a target="_blank" href="https://doi.org/10.1016/j.resmic.2022.103936">10.1016/j.resmic.2022.103936</a>. Epub ahead of print).</p>
