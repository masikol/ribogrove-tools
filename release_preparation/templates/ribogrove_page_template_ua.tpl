<!DOCTYPE html>
<html lang="uk">
<head>
  <title>Дзеркало RiboGrove</title>
  <link rel="stylesheet" href="css/style.css">
  <link rel="stylesheet" href="css/release-page.css">
</head>

<body>
<div class="container">
<div class="lang-section">
<p>&#127760;<br><a href="index.html">English</a> / <a href="ribogrove_{{ ribogrove_release_number }}_be.html">Беларуская</a> / <span class="active-lang">Українська</span> / <a href="ribogrove_{{ ribogrove_release_number }}_ru.html">Русский</a></p>
</div>
<p><div class="warn-section"><a href="http://mbio.bas-net.by/cager/en/ribogrove">Головний сайт</a>, де ми розміщуємо базу даних RiboGrove, може бути недосяжним поза межами Білорусі через технічні проблеми та загальну біду. Актуальний реліз RiboGrove можна знайти тут &mdash; нехай база даних буде досяжна принаймни такою мірою.</div></p><br>
<h1 style="display:none">RiboGrove</h1>
<img src="img/RiboGrove_logo.png"><br><br>
<h2>Зміст</h2>
<ul>
<li><a href="#overview">Що таке RiboGrove</a>
<ul style="list-style-type: none;">
<li><a href="#other-databases">RiboGrove і інші бази даних 16S рРНК</a></li>
<li><a href="#genome-categories">Категорії геномів</a></li>
</ul>
</li>
<li><a href="#citing-ribogrove">Цитування RiboGrove</a></li>
<li><a href="#downloads">Файли</a>
<ul style="list-style-type: none;">
  <li><a href="#latest-release">Актуальний реліз RiboGrove</li>
  <li><a href="#release-archive">Архів релізів RiboGrove</li>
</ul>
</li>
<li><a href="#summary">Статистичний опис</a>
<ul style="list-style-type: none;">
<li><a href="#db-size">Розмір RiboGrove</a></li>
<li><a href="#gene-lens">Довжина генів 16S рРНК</a></li>
<li><a href="#copy-number">Кількість копій генів 16S рРНК</a></li>
<li><a href="#top-longest">Топ-10 найдовших генів 16S рРНК</a></li>
<li><a href="#top-shortest">Топ-10 найкоротших генів 16S рРНК</a></li>
<li><a href="#top-copy-num">Топ-10 геномів за найбільшею кілкістю копій генів 16S рРНК</a></li>
<li><a href="#top-var">Топ-10 геномів за найбільшею мінливістю генів 16S рРНК</a></li>
</ul>
</li>
<li><a href="#searching-data">Пошук даних у RiboGrove</a>
<ul style="list-style-type: none;">
  <li><a href="#header-format">Формат заголовку</a></li>
<li><a href="#subset-by-header">Вибір послідовностей по заголовку</a></li>
<li><a href="#subset-by-len">Вибір послідовностей по довжині</a></li>
<li><a href="#select-head">Вибір даних заголовків</a></li>
</ul>
</li>
</ul>
<hr>
<div id="overview" class="pad-anchor"></div>
<h2>Що таке RiboGrove</h2>
<p>RiboGrove &mdash; це база даних послідовностей генів 16S рРНК бактерій та архей.</p>
<p>RiboGrove базується на базі даних <a href="https://ncbi.nlm.nih.gov/refseq/">RefSeq</a>. RiboGrove утримує тільки повнорозмірні послідовності генів 16S рРНК, а самі послідовності походять із повністю зібраних геномів прокаріотів, депонуваних у RefSeq. Тому ми постулюємо високу надійність послідовностей, розміщених у RiboGrove.</p>
<div id="other-databases" class="pad-anchor"></div>
<h3>RiboGrove і інші бази даних 16S рРНК</h3>
<p>Таблиця нижче показує, чим (якісно) RiboGrove відрізнювається від подібних баз даних послідовностей РНК, а саме <a href="https://www.arb-silva.de/">Silva</a>, <a href="http://rdp.cme.msu.edu/">RDP</a> та <a href="https://greengenes.secondgenome.com/">Greengenes</a>. Кажучи коротко, RiboGrove менша за аналоги щодо кількості послідовностей та їх різноманітності, але виграє у надійності послідовностей.</p>
<table class="sum-table">
<tbody>
  <tr class="sumtab-row"><th></th><th class="alcentr">RiboGrove</th><th class="alcentr">Silva</th><th class="alcentr">RDP</th><th class="alcentr">Greengenes</th></tr>
  <tr class="sumtab-row"><td class="feature-col">Репрезентовані організми</td><td class="alcentr">Бактерії<br>Археї</td><td class="alcentr">Бактерії<br>Археї<br>Еукаріоти</td><td class="alcentr">Бактерії<br>Археї<br>Еукаріоти</td><td class="alcentr">Бактерії<br>Археї</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Репрезентовані субодиниці рибосом</td><td class="alcentr">Мала</td><td class="alcentr">Велика<br>Мала</td><td class="alcentr">Велика<br>Мала</td><td class="alcentr">Мала</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Утримує послідовності зібраних геномів</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Утримує послідовності ампліконів</td><td class="alcentr">Ні</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Утримує неповні послідовності генів</td><td class="alcentr">Ні</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
</tbody>
</table>
<div id="genome-categories" class="pad-anchor"></div>
<h3>Категорії геномів</h3>
<p>Всі геноми, з якіх бралися дані для створення RiboGrove, були розділені на три категорії згідно з їхньою надійностю:</p>
<ol style="list-style:none;">
  <li><strong>Категорія 1</strong> (найбільша надійність). Геноми, які не демонструють прикмет низкоякісної збірки і які були секвенувані за допомогою технології PacBio або комбінації (Oxford&nbsp;Nanopore&nbsp;+&nbsp;Illumina).</li>
  <li><strong>Категорія 2.</strong> Геноми, які не демонструють прикмет низкоякісної збірки, секвенувані за допомогою будь-якої іншої технології (або ті, для яких технологія секвенування не була задекларувана).</li>
  <li><strong>Категорія 3</strong> (найменша надійність). Геноми, які демонструють принаймні одну прикмету низкоякісної збірки.</li>
</ol>
<p>Прикмети низкоякісної збірки</p>
<ul>
  <li>Геном утримує принаймні одну звироднілу основу у послідовностях його генів 16S рРНК.</li>
  <li>Збірка генома утримує принаймні один запис бази даних RefSeq, у назві якої є фраза <q>map unlocalized</q>, і цей запис утримує послідовность гена 16S рРНК або його частину.</li>
</ul>
<br>
<p>Програми, за допомогую яких RiboGrove була створена, находится у наступному ҐітХаб-репозиторії: <a href="https://github.com/masikol/ribogrove-tools">ribogrove-tools</a>.</p>
<hr>
<div id="citing-ribogrove" class="pad-anchor"></div>
<h2>Цитування RiboGrove</h2>
<p>Якщо база даних RiboGrove сталася корисною для ваших досліджень, посилайтеся, будь ласка, на наступну публікацію:</p>
<p>Maxim A. Sikolenko, Leonid N. Valentovich. &ldquo;RiboGrove: a database of full-length prokaryotic 16S rRNA genes derived from completely assembled genomes&rdquo; // Research in Microbiology, Volume 173, Issue 4, May 2022, 103936.<br> (ЦІО: <a href="https://doi.org/10.1016/j.resmic.2022.103936">10.1016/j.resmic.2022.103936</a>).</p>
<hr>
<div id="downloads" class="pad-anchor"></div>
<h2>Завантаження</h2>
<div id="latest-release" class="pad-anchor"></div>
<h3>Актуальний реліз RiboGrove &mdash; {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
<p>Реліз заснований на базі даних RefSeq {{ refseq_release }}.</p>
<ul>
<li>Fasta-файл повнорозмірних послідовностей генів 16S рРНК. <a href="ribogrove_releases/{{ ribogrove_release_number }}/ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Завантажити (&bdquo;ґзіп&rsquo;нутий&ldquo; fasta-файл, {{ final_fasta_fsize_fmt }} Мб)</a></li>
<li>&bdquo;Сира&ldquo; версія RiboGrove &mdash; у неї часткові послідовності генів залишені. <a href="ribogrove_releases/{{ ribogrove_release_number }}/raw_ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Завантажити (&bdquo;ґзіп&rsquo;нутий&ldquo; fasta-файл, {{ raw_fasta_fsize_fmt }} Мб)</a></li>
<li>Метадані. <a href="ribogrove_releases/{{ ribogrove_release_number }}/metadata_ribogrove_{{ ribogrove_release_number }}.zip">Завантажити (zip-архів {{ metadata_fsize_fmt }} Мб)</a><br>
<details>
  <summary><small>Яка саме інформація утримується у метаданих?</small></summary>
  <small>
    <p>Метадані складаються з наступних файлів:</p>
    <ol>
    <li>
      <span class="samp-highl">source_RefSeq_genomes.tsv</span><br>
      Це TSV-файл, який утримує інформацію про те, які геноми выкористовувалися для створення RiboGrove.
    </li>
    <li>
      <span class="samp-highl">gene_seqs_statistics.tsv</span>, <span class="samp-highl">raw_gene_seqs_statistics.tsv</span><br>
      Це TSV-файли, які утримують інформацію пра нуклеотидний склад, розмір, геномну і таксономічну приналежність послідовностей генів. &bdquo;Сира&ldquo; (&bdquo;raw&ldquo;) версія додатково утримує інформацію про неповні гени.
    </li>
    <li>
      <span class="samp-highl">per_replicon_statistics.tsv</span>, <span class="samp-highl">raw_per_replicon_statistics.tsv</span><br>
      Це TSV-файли, які утримують інформацію про кількісць генів 16S рРНК у кожної геномної послідовності з БД RefSeq і про розміри цих генів. &bdquo;Сира&ldquo; (&bdquo;raw&ldquo;) версія додатково утримує інформацію про неповні гени.
    </li>
    <li>
      <span class="samp-highl">categories.tsv</span><br>
      Це TSV-файл, який утримує інформацію про категорію кожного геному, який був використаний для створення RiboGrove. Також у цьому файлі є інформація про технології, які выкористовувалися для секвенування цих геномів.
    </li>
    <li>
      <span class="samp-highl">taxonomy.tsv</span><br>
      Це TSV-файл, який утримує таксономічну приналежність кожного геному і гена.
    </li>
    <li>
      <span class="samp-highl">intragenic_repeats.tsv</span><br>
      Це TSV-файл, який утримує інформацію про внутрішьогенні повтори у послідовностях генів, знайдені за допомогою програми <a href="https://github.com/deprekate/RepeatFinder">RepeatFinder</a>.
    </li>
    <li>
      <span class="samp-highl">cmscan_output_table.tblout</span><br>
      Це TSV-файл, який утримує повний результат роботи програми <a href="http://eddylab.org/infernal/">cmscan</a>, після оброблення всіх послідовностей генів 16S рРНК.
    </li>
    <li>
      <span class="samp-highl">entropy_summary.tsv</span><br>
      Це TSV-файл, який утримує інформацію про внутрішньогеномну мінливість генів 16S рРНК. Внутрішньогеномну мінливість розраховували тільки для геномів категорії&nbsp;1, у яких є принаймні два гени 16S рРНК. Для оцінки мінливості выкористовували ентропію Шенона: вирівнювали послідовності генів кожного геному за допомогою програми <a href="https://drive5.com/muscle5/">MUSCLE</a>, і після розраховували ентропію для кожної колонки множинного вирівнювання.
    </li>
    <li>
      <span class="samp-highl">QIIME2-compatible-taxonomy.txt</span><br>
      Це TSV-файл, який можна використовувати для того, щоб навчити класіфікатор QIIME2 (гл. <a href="https://docs.qiime2.org">інструкцію</a> ).
    </li>
  </ol></small>
</details>
</li>
</ul>
<p>Файли fasta стиснуті за допомогою програми gzip, а файли метаданих стиснуті у zip-архів. Щоби розпакувати їх, користувачі Linux та Mac OS можуть скористатися саме програмамі gzip та zip; оці програми звичайно установлені разом із операційною системою. Користувачам Windows може придатися безкоштовна програма-архіватор <a href="https://www.7-zip.org/">7-Zip</a>.</p>
<div id="release-archive" class="pad-anchor"></div>
<h3>Архів релізів RiboGrove</h3>
<p>Всі попередні релізи можна знайти в <a href="archive/ribogrove_archive_en.html"> архіві релізів RiboGrove</a>. На жаль, він ведеться тільки англійською.</p>
<hr>
<div id="summary" class="pad-anchor"></div>
<h2>Статистичний опис</h2>
<div id="db-size" class="pad-anchor"></div>
<table class="sum-table">
<caption>Розмір RiboGrove</caption>
<tbody>
<tr><th></th><th class="numcol">Бактерії</th><th class="numcol">Археї</th><th class="numcol">Усього</th></tr>
<tr class="sumtab-row">
<td>Кількість послідовностей генів</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Кількість унікальних послідовностей</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Кількість видів</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Кількість геномів</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Кількість геномів 1-ої категорії</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Кількість геномів 2-ої категорії</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Кількість геномів 3-ої категорії</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Total'] }}</td>
</tr>
</tbody>
</table>
<div id="gene-lens" class="pad-anchor"></div>
<table class="sum-table">
<caption>Довжина генів 16S рРНК</caption>
<tbody>
<tr><th></th><th class="numcol">Бактерії</th><th class="numcol">Археї</th></tr>
<tr class="sumtab-row">
<td>Мінімальна довжина, п.о.</td>
<td class="numcol">{{ ribogrove_len_dict['min']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['min']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>25-й перцентиль, п.о. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['25perc']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['25perc']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Медіана, п.о. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['median']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['median']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>75-й перцентиль, п.о. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['75perc']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['75perc']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Середня довжина, п.о. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['mean']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['mean']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Мода, п.о. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['modes']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['modes']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Максимальна довжина, п.о.</td>
<td class="numcol">{{ ribogrove_len_dict['max']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['max']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Сереньоквадратичне відхилення, п.о.<sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['std']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['std']['Archaea'] }}</td>
</tr>
</tbody>
</table>
<p><sup>*</sup> Метрики, позначені зірочкою, розрахувано зі попередньою нормалізацією: вирахуванням <em>медіан</em> довжин генів у межах виду.</p>
<div id="copy-number" class="pad-anchor"></div>
<table class="sum-table">
<caption>Кількість копій генів 16S рРНК</caption>
<tbody>
<tr><th class="numcol">Кількість копій <sup>*</sup></th><th class="numcol">Кількість видів</th><th class="numcol">Відсоток видів, %</th></tr>
{% for _, row in ribogrove_copy_number_df.iterrows() %}
<tr class="sumtab-row">
<td class="numcol">{{ row['copy_number'] }}</td>
<td class="numcol">{{ row['number_of_species'] }}</td>
<td class="numcol">{{ row['percent_of_species'] }}</td>
</tr>
{% endfor %}
</tbody>
</table>
<p><sup>*</sup> Медіанна кількість копій у межах виду.</p>

<div id="top-longest" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 найдовших генів 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Організм</th><th class="numcol">Довжина, п.о.</th><th class="alnleft">SeqID <sup>*</sup></th><th class="alnleft">Ідентифікатор<br> геномної збірки</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерії</td>
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
<td colspan="4" class="subhead">Археї</td>
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
<p><sup>*</sup> SeqID &mdash; це ідентифікатор послідовності RiboGrove.</p>

<div id="top-shortest" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 найкоротших генів 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Організм</th><th class="numcol">Довжина, п.о.</th><th class="alnleft">SeqID <sup>*</sup></th><th class="alnleft">Ідентифікатор<br> геномної збірки</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерії</td>
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
<td colspan="4" class="subhead">Археї</td>
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
<p><sup>*</sup> SeqID &mdash; це ідентифікатор послідовності RiboGrove.</p>

<div id="top-copy-num" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 геномів за найбільшею кілкістю копій генів 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Організм</th><th class="numcol">Кількість копій</th><th class="alnleft">Ідентифікатор<br> геномної збірки</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерії</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['title']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Археї</td>
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
<table class="sum-table">
<caption>Топ-10 геномів за найбільшею мінливістю генів 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Організм</th><th class="numcol">Сума ентропії<sup> *</sup>, біти</th><th class="numcol">Середня ентропія<sup> *</sup>, біти</th><th class="numcol">Кількість мінливих позіцій</th><th class="numcol">Кількість копій</th><th class="alnleft">Ідентифікатор<br> геномної збірки</th></tr>
<tr>
<td colspan="6" class="subhead">Бактерії</td>
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
<td colspan="6" class="subhead">Археї</td>
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

<p><sup>*</sup> Ентропія тут &mdash; це <a href="https://uk.wikipedia.org/wiki/Інформаційна_ентропія">ентропія Шеннона</a>, яку розрахувовали для кожної колонки множинного вирівнювання повнорозмірних послідовностей генів 16S рРНК кожного геному. І далі ентропію підсумовували (колонка <q>Сума ентропії</q>) і усереднювали (колонка <q>Середня ентропія</q>).</p>
<p><sup>**</sup> <i>Halomicrobium</i> sp. ZPS1 &mdash; показовий випадок. У цьому геномі два гени 16S рРНК, тому ентропія є рівною кількості нуклеотидів, які не співпадають у їхніх послідовностях. Коли підрахувати відсоток ідентичності (скорочено ВІ) послідовностей оцих двох генів, отримуємо 90,70%! Нагадуємо, що за поріг ВІ послідовностей генів 16S рРНК для віднесення організмів до відрізных <em>родів</em> звычайно береться (хоча і довільно) значення 95%.</p>
<hr>
<div id="searching-data" class="pad-anchor"></div>
<h2>Пошук даних у RiboGrove</h2>
<p>RiboGrove &mdash; це дуже мінімалістична база даних: вона складається з двох файлів fasta із додатковими метаданими. Отже, розширених засобів пошуку немає. Ми визнаємо цю проблему і залишаємо нижче делілька рекомендацій щодо пошуку даних у RiboGrove.</p>
<div id="header-format" class="pad-anchor"></div>
<h3>Формат заголовку</h3>
<p>Записи послідовностей RiboGrove у форматі fasta мають наступний формат заголовку:</p>
<p class="samp-highl samp-vwide">&gt;NZ_CP079719.1:86193-87742_plus Bacillus_velezensis ;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus; category:2</p>
<p>Основні частини заголовку відокремлені один від одного пробілами. В заголовку повинні бути чотири такі частини:</p>
<ol>
<li>Ідентифікатор послідовності (SeqID): <span class="samp-highl">NZ_CP079719.1:86193-87742_plus</span>. SeqID, у свою чергу, складається з трох частинок:<ol>
<li>код доступу послідовності RefSeq, від якої походить послідовність гену: <span class="samp-highl">NZ_CP079719.1</span>;</li>
<li>координати гену у цієї геномної послідовності RefSeq: <span class="samp-highl">86193-87742</span> (нумерація йде з 1-ці, інтервали включають ліві та праві межі);</li>
<li>ланцюг геномної послідовності, на якій ген знаходиться: <span class="samp-highl">plus</span> (або <span class="samp-highl">minus</span>).</li>
</ol></li>
<li>Повне імя організму, згідно з <a href="https://ncbi.nlm.nih.gov/taxonomy">базою даних NCBI Taxonomy</a>: <span class="samp-highl">Bacillus_velezensis</span>.</li>
<li>Рядок таксономії, а саме назви домену (<span class="samp-highl">Bacteria</span>), типу (<span class="samp-highl">Firmicutes</span>), класу (<span class="samp-highl">Bacilli</span>), ряду (<span class="samp-highl">Bacillales</span>), родини (<span class="samp-highl">Bacillaceae</span>) та ріду (<span class="samp-highl">Bacillus</span>). Назви відокремлені та оточені крапками з комою (<span class="samp-highl">;</span>).</li>
<li>Категорія геному, із якого послідовність гену походить: (<span class="samp-highl">category:2</span>). Геноми категорії 1 є найнадійними, а категорії 3 &mdash; найменш надійними.</li>
</ol>
<h3>Вибір послідовностей</h3>
<p>Відбирати ті чи інші послідовності із файлів fasta можна за допомогою програми <strong>Seqkit</strong> (<a href="https://github.com/shenwei356/seqkit">ҐітХаб-репозиторій</a>, <a href="https://bioinf.shenwei.me/seqkit/">документація</a>). Це безкоштовна, кросплатформенна, багатофункціональна та доволі швидка програма, яка може обробляти як &bdquo;ґзіп&rsquo;нуті&ldquo;, так і прості нестиснуті файли fasta. Для відбіру послідовностей зручно выкористовувати програми <span class="samp">seqkit grep</span> та <span class="samp">seqkit seq</span>.</p>
<div id="subset-by-header" class="pad-anchor"></div>
<h4>Вибір послідовностей по заголовку</h4>
<p>Нехай ми робимо вибірку з файла <span class="samp">ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</span>. Корисними можуть статися наступні приклади команд користання програмою <span class="samp">seqkit grep</span>:</p>
<p><strong>Приклад 1</strong>. Вибрати одну послідовність по SeqID.</p>
<p class="samp-highl samp-vwide">seqkit grep -p "NZ_CP079719.1:86193-87742_plus" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опція <span class="samp">-p</span> задає рядок, який програма буде шукати у fasta-заголовках (на самій справі, тільки в тій частині заголовків, де запісаний SeqID).</p>
<p><strong>Приклад 2</strong>. Вибрати всі послідовності генів однієї геномної послідовності RefSeq по коду доступу.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "NZ_CP079719.1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Тут потрібні ще дві опціі: <span class="samp">-n</span> та <span class="samp">-r</span>. Перша вказує програмі шукати підрядки в заголовках цілком, а не тільки в ідентифікаторах (SeqID). Друга опція вказує програмі, що їй треба шукати не тільки ті заголовки, які цілком співпадають із шуканим рядком, а також і ті, які утримують шуканий рядок як свою частину.</p>
<p><strong>Приклад 3</strong>. Вибрати всі послідовності генів одного геному (ідентифікатор геномної збірки <a href="https://ncbi.nlm.nih.gov/assembly/10577151">10577151</a>), який має два реплікони: <a href="https://ncbi.nlm.nih.gov/nuccore/NZ_CP079110.1">NZ_CP079110.1</a> та <a href="https://ncbi.nlm.nih.gov/nuccore/NZ_CP079111.1">NZ_CP079111.1</a>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nr -p "NZ_CP079110.1" -p "NZ_CP079111.1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Приклад 4</strong>. Вибрати всі послідовності актинобактерий.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ";Actinobacteria;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Про всяк випадок, оточюйте таксономію крапками з комою (<span class="samp-highl">;</span>).</p>
<p><strong>Приклад 5</strong>. Вибрати всі послідовності з геномів категорії 1.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "category:1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Приклад 6</strong>. Вибрати всі послідовності крім тих, які належать до типу <i>Firmicutes</i>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nvrp ";Firmicutes;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Зверніть увагу на опцію <span class="samp">-v</span> у рядку з опціями <span class="samp">-nvrp</span>. Оця опція інвертує поведінку програми, і тепер вона знайде всі послідовності, заголовки яких <em>не</em> утримують рядку <span class="samp">;Firmicutes;</span>.</p>
<div id="subset-by-len" class="pad-anchor"></div>
<h4>Вибір послідовностей по довжині</h4>
<p>Відбірати послідовності по довжині можна за допомогою програми <span class="samp">seqkit seq</span>.</p>
<p><strong>Приклад 1</strong>. Вибрати всі послідовності довше ніж 1600 п.о.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1601 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опція <span class="samp">-m</span> задає мінімальну довжину послідовностей, які програма подасть на вихід.</p>
<p><strong>Приклад 2</strong>. Вибрати всі послідовності коротше ніж 1500 п.о.</p>
<p class="samp-highl samp-vwide">seqkit seq -M 1499 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опція <span class="samp">-M</span> задає максимальну довжину послідовностей, які програма подасть на вихід.</p>
<p><strong>Приклад 3</strong>. Вибрати всі послідовності, довжина яких знаходится в інтервалі [1500, 1600] п.о.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1500 -M 1600 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<div id="select-head" class="pad-anchor"></div>
<h4>Вибір даних заголовків</h4>
<p>Інколи буває корисно одержати тільки дані заголовків із файла fasta. Це можна робити за допомогою програмы <span class="samp">seqkit seq</span>.</p>
<p><strong>Приклад 1</strong>. Вибрати просто всі заголовки.</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опція <span class="samp">-n</span> каже програмі, щоби вона подавала на вихід тільки заголовки.</p>
<p><strong>Приклад 2</strong>. Вибрати всі ідентифікатори послідовностей (тобто частини заголовків до першого пробіла).</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опція <span class="samp">-i</span> каже програмі, щоби вона подавала на вихід тільки ідентифікатори.</p>
<p><strong>Приклад 3</strong>. Вибрати всі (RefSeq&rsquo;овські) коди доступу.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f1 -d':' | sort | uniq</p>
<p class="samp-comment">Така команда спрацює, тільки якщо на комп&rsquo;ютері установлені утиліти <span class="samp">sort</span>, <span class="samp">cut</span> та <span class="samp">uniq</span> (на Linux та Mac OS оці програми звичайно встановлені разом із сістемою).</p>
<br>
<hr>
<p>RiboGrove, {{ ribogrove_release_date }}</p>
</div>
</body>
</html>
