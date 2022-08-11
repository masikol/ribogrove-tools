<!DOCTYPE html>
<html lang="be">
<head>
  <title>Люстэрка RiboGrove</title>
  <link rel="stylesheet" href="css/style.css">
  <link rel="stylesheet" href="css/release-page.css">
</head>

<body>
<div class="container">
<div class="lang-section">
<p>&#127760;<br><a href="index.html">English</a> / <span class="active-lang">Беларуская</span> / <a href="ribogrove_{{ ribogrove_release_number }}_ua.html">Українська</a> / <a href="ribogrove_{{ ribogrove_release_number }}_ru.html">Русский</a></p>
</div>
<p><div class="warn-section"><a href="http://mbio.bas-net.by/cager/en/ribogrove">Асноўны сайт</a>, дзе размяшчаецца база дадзеных RiboGrove, можа быць недасяжным па-за межамі Беларусі праз тэхнічныя прычыны ды агульную катастрофу. Таму існуе гэтае люстэрка, а файлы RiboGrove дасяжныя праз Дропбокс (спасылкі <a href="#downloads">ніжэй</a>).</div></p><br>
<h1 style="display:none">RiboGrove</h1>
<img src="img/RiboGrove_logo.png"><br><br>
<h2>Змест</h2>
<ul>
<li><a href="#overview">Што такое RiboGrove</a>
<ul class="ribogrove-nested-list">
<li><a href="#other-databases">RiboGrove і іншыя базы дадзеных 16S рРНК</a></li>
<li><a href="#genome-categories">Катэгорыі геномаў</a></li>
</ul>
</li>
<li><a href="#citing-ribogrove">Цытаванне RiboGrove</a></li>
<li><a href="#downloads">Файлы</a>
<ul class="ribogrove-nested-list">
<li><a href="#current-release">Актуальны выпуск RiboGrove</li>
<li><a href="#release-archive">Архіў выпускаў RiboGrove</li>
</ul>
</li>
<li><a href="#summary">Статыстычнае апісанне</a>
<ul class="ribogrove-nested-list">
<li><a href="#db-size">Памер RiboGrove</a></li>
<li><a href="#gene-lens">Даўжыня генаў 16S рРНК</a></li>
<li><a href="#copy-number">Колькасть копій генаў 16S рРНК</a></li>
<li><a href="#top-longest">Топ-10 найдаўжэйшых генаў 16S рРНК</a></li>
<li><a href="#top-shortest">Топ-10 найкарацейшых генаў 16S рРНК</a></li>
<li><a href="#top-copy-num">Топ-10 геномаў з найбольшай колькасцю копій генаў 16S рРНК</a></li>
<li><a href="#top-var">Топ-10 геномаў з найбольшай зменлівасцю генаў 16S рРНК</a></li>
</ul>
</li>
<li><a href="#searching-data">Пошук дадзеных у RiboGrove</a>
<ul class="ribogrove-nested-list">
<li><a href="#header-format">Фармат загалоўка</a></li>
<li><a href="#subset-by-header">Адбор паслядоўнасцяў паводле загалоўкаў</a></li>
<li><a href="#subset-by-len">Адбор паслядоўнасцяў паводле даўжыні</a></li>
<li><a href="#select-head">Адбор дадзеных загалоўкаў</a></li>
</ul>
</li>
</ul>
<hr>
<div id="overview" class="pad-anchor"></div>
<h2>Што такое RiboGrove</h2>
<p>RiboGrove &mdash; гэта база дадзеных паслядоўнасцей генаў 16S рРНК бактэрый і архей.</p>
<p>RiboGrove асноўваецца на базе дадзеных <a href="https://ncbi.nlm.nih.gov/refseq/">RefSeq</a>. RiboGrove змяшчае толькі поўнапамерныя паслядоўнасці генаў 16S рРНК, а самі паслядоўнасці паходзяць з цалкам сабраных геномаў прокарыёт, дэпанаваных у RefSeq. Таму мы пастулюем высокую надзейнасць паслядоўнасцей, размешчаных у RiboGrove.</p>
<div id="other-databases" class="pad-anchor"></div>
<h3>RiboGrove і іншыя базы дадзеных 16S рРНК</h3>
<p>Табліца ніжэй паказвае, чым (якасна) RiboGrove адрозніваецца ад падобных баз дадзеных паслядоўнасцей РНК, а менавіта <a href="https://www.arb-silva.de/">Silva</a>, <a href="http://rdp.cme.msu.edu/">RDP</a> і <a href="https://greengenes.secondgenome.com/">Greengenes</a>. Кажучы каротка, RiboGrove саступае аналагам у колькасці паслядоўнасцей і іх разнастайнасці, але пераўзыходзіць у надзейнасці паслядоўнасцей.</p>
<table class="sum-table">
<tbody>
  <tr class="sumtab-row"><th></th><th class="alcentr">RiboGrove</th><th class="alcentr">Silva</th><th class="alcentr">RDP</th><th class="alcentr">Greengenes</th></tr>
  <tr class="sumtab-row"><td class="feature-col">Прадстаўленыя арганізмы</td><td class="alcentr">Бактэрыі<br>Археі</td><td class="alcentr">Бактэрыі<br>Археі<br>Эўкарыёты</td><td class="alcentr">Бактэрыі<br>Археі<br>Эўкарыёты</td><td class="alcentr">Бактэрыі<br>Археі</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Прадстаўленыя субадзінкі рыбасом</td><td class="alcentr">Малая</td><td class="alcentr">Вялікая<br>Малая</td><td class="alcentr">Вялікая<br>Малая</td><td class="alcentr">Малая</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Змяшчае паслядоўнасці сабраных геномаў</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Змяшчае паслядоўнасці ампліконаў</td><td class="alcentr">Не</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Змяшчае няпоўныя паслядоўнасці генаў</td><td class="alcentr">Не</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
</tbody>
</table>
<div id="genome-categories" class="pad-anchor"></div>
<h3>Катэгорыі геномаў</h3>
<p>Усе геномы, з якіх браліся дадзеныя для стварэння RiboGrove, былі падзеленыя на тры катэгорыі паводле іхняй надзейнасці:</p>
<ol style="list-style:none;">
  <li><strong>Катэгорыя 1</strong> (найбольшая надзейнасць). Геномы, якія не дэманструюць прыкмет нізкаякаснай зборкі і якія былі секвенаваныя з дапамогай тэхналогіі PacBio альбо камбінацыі (Oxford&nbsp;Nanopore&nbsp;+&nbsp;Illumina).</li>
  <li><strong>Катэгорыя 2.</strong> Геномы, якія не дэманструюць прыкмет нізкаякаснай зборкі, секвенаваныя з дапамогай любой іншай тэхналогіі (або тыя, для якіх тэхналогія секвенавання не была ўказаная).</li>
  <li><strong>Катэгорыя 3</strong> (найменшая надзейнасць). Геномы, якія дэманструюць хаця б адну прыкмету нізкаякаснай зборкі.</li>
</ol>
<p>Прыкметы нізкаякаснай зборкі:</p>
<ul>
  <li>Геном утрымлівае хаця б адну выраджаную аснову у паслядоўнасцях яго генаў 16S рРНК.</li>
  <li>Зборка генома утрымлівае хаця б адзін запіс базы дадзеных RefSeq, у назве якой ёсць фраза <q>map unlocalized</q>, і гэты запіс утрымлівае паслядоўнасць гена 16S рРНК або яго частку.</li>
</ul>
<br>
<p>Праграмы, з дапамогай якіх была створана RiboGrove, знаходзяцца ў наступным ГітХаб-рэпазіторыі: <a href="https://github.com/masikol/ribogrove-tools">ribogrove-tools</a>.</p>
<hr>
<div id="citing-ribogrove" class="pad-anchor"></div>
<h2>Цытаванне RiboGrove</h2>
<p>Калі база дадзеных RiboGrove прыдалася вам у вашых даследаваннях, спасылайцеся, калі ласка, на наступную публікацыю:</p>
<p>Maxim A. Sikolenko, Leonid N. Valentovich. &ldquo;RiboGrove: a database of full-length prokaryotic 16S rRNA genes derived from completely assembled genomes&rdquo; // Research in Microbiology, Volume 173, Issue 4, May 2022, 103936.<br> (<span style="cursor:help"title="Лічбавы ідэнтыфікатар аб'екта">ЛІА</span>: <a href="https://doi.org/10.1016/j.resmic.2022.103936">10.1016/j.resmic.2022.103936</a>).</p>
<hr>
<div id="downloads" class="pad-anchor"></div>
<h2>Файлы</h2>
<div id="latest-release" class="pad-anchor"></div>
<h3>Актуальны выпуск RiboGrove &mdash; {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
<p>Выпуск заснованы на базе дадзеных RefSeq {{ refseq_release }}.</p>
<ul>
<li>Fasta-файл з поўнапамернымі паслядоўнасцямі генаў 16S рРНК. <a href="ribogrove_releases/{{ ribogrove_release_number }}/ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Спампаваць (&bdquo;гзіп&rsquo;нуты&ldquo; fasta-файл, {{ final_fasta_fsize_fmt }} Мб)</a></li>
<li>&bdquo;Сырая&ldquo; версія RiboGrove &mdash; яна змяшчае ў тым ліку і частковыя паслядоўнасці генаў. <a href="ribogrove_releases/{{ ribogrove_release_number }}/raw_ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Спампаваць (&bdquo;гзіп&rsquo;нуты&ldquo; fasta-файл, {{ raw_fasta_fsize_fmt }} Мб)</a></li>
<li>Метададзеныя. <a href="ribogrove_releases/{{ ribogrove_release_number }}/metadata_ribogrove_{{ ribogrove_release_number }}.zip">Спампаваць (zip-архіў {{ metadata_fsize_fmt }} Мб)</a><br>
<details>
  <summary><small>Якая менавіта інфармацыя утрымліваецца ў метададзеных?</small></summary>
  <small>
    <p>Метададзеныя складаюцца з наступных файлаў:</p>
    <ol>
    <li>
      <span class="samp-highl">source_RefSeq_genomes.tsv</span><br>
      Гэта TSV-файл, які змяшчае інфармацыю пра тое, якія геномы выкарыстоўваліся для стварэння RiboGrove.
    </li>
    <li>
      <span class="samp-highl">gene_seqs_statistics.tsv</span>, <span class="samp-highl">raw_gene_seqs_statistics.tsv</span><br>
      Гэта TSV-файлы, якія змяшчаюць інфармацыю пра нуклеатыдны склад, памер, геномную і таксанамічную прыналежнасць паслядоўнасцей генаў. &bdquo;Сырая&ldquo; (&bdquo;raw&ldquo;) версія дадаткова ўтрымлівае інфармацыю пра няпоўныя гены.
    </li>
    <li>
      <span class="samp-highl">per_replicon_statistics.tsv</span>, <span class="samp-highl">raw_per_replicon_statistics.tsv</span><br>
      Гэта TSV-файлы, якія змяшчаюць інфармацыю пра колькасць генаў 16S рРНК у кожнай геномнай паслядоўнасці з БД RefSeq і пра памеры гэтых генаў. &bdquo;Сырая&ldquo; (&bdquo;raw&ldquo;) версія дадаткова ўтрымлівае інфармацыю пра няпоўныя гены.
    </li>
    <li>
      <span class="samp-highl">categories.tsv</span><br>
      Гэта TSV-файл, які змяшчае інфармацыю пра катэгорыю кожнага геному, які быў выкарыстаны для стварэння RiboGrove. Таксама у гэтым файле ёсць інфамацыя пра тэхналогіі, якія выкарыстоўваліся для секвенавання гэтых геномаў.
    </li>
    <li>
      <span class="samp-highl">taxonomy.tsv</span><br>
      Гэта TSV-файл, які змяшчае таксанамічную прыналежнасць кожнага геному і гена.
    </li>
    <li>
      <span class="samp-highl">intragenic_repeats.tsv</span><br>
      Гэта TSV-файл, які змяшчае інфармацыю пра унутрыгенныя паўторы ў паслядоўнасцях генаў, знойдзеныя з дапамогай праграмы <a href="https://github.com/deprekate/RepeatFinder">RepeatFinder</a>.
    </li>
    <li>
      <span class="samp-highl">cmscan_output_table.tblout</span><br>
      Гэта TSV-файл, які змяшчае поўны вынік работы праграмы <a href="http://eddylab.org/infernal/">cmscan</a>, пасля апрацоўкі ўсіх паслядоўнасцяў генаў 16S рРНК.
    </li>
    <li>
      <span class="samp-highl">entropy_summary.tsv</span><br>
      Гэта TSV-файл, які змяшчае інфармацыю пра унутрыгеномную зменлівасць генаў 16S рРНК. Унутрыгеномную зменлівасць вылічалі толькі для геномаў катэгорыі&nbsp;1, у якіх ёсць хаця б два гены 16S рРНК. Для ацэнкі зменлівасці выкарыстоўвалі энтрапію Шэнана: выраўноўвалі паслядоўнасці генаў кожнага геному з дапамогай праграмы <a href="https://drive5.com/muscle5/">MUSCLE</a>, і пасля вылічалі энтрапію для кожнай калонкі множнага выраўноўвання.
    </li>
    <li>
      <span class="samp-highl">QIIME2-compatible-taxonomy.txt</span><br>
      Гэта TSV-файл, які можна выкарыстоўваць для таго, каб навучыць класіфікатар QIIME2 (гл. <a href="https://docs.qiime2.org">інструкцыю</a> ).
    </li>
  </ol></small>
</details>
</li>
</ul>
<p>Файлы fasta сціснутыя з дапамогаю праграмы gzip, а файлы з метададзенымі сціснутыя ў zip-архіў. Каб распакаваць іх, карыстальнікі Linux і Mac OS могуць скарастыцца якраз праграмамі gzip і zip; гэтыя праграмы звычайна ўсталяваныя разам з аперыцыйнай сістэмай. Карыстальнікам Windows можа прыдацца бясплатная праграма-архіватар <a href="https://www.7-zip.org/">7-Zip</a>.</p>
<div id="release-archive" class="pad-anchor"></div>
<h3>Архіў выпускаў RiboGrove</h3>
<p>Усе папярэднія выпускі можна знайсці ў <a href="archive/ribogrove_archive_en.html"> архіве выпускаў RiboGrove</a>. На жаль, ён вядзецца толькі па-англійску.</p>
<hr>
<div id="summary" class="pad-anchor"></div>
<h2>Статыстычнае апісанне</h2>
<div id="db-size" class="pad-anchor"></div>
<table class="sum-table">
<caption>Памер RiboGrove</caption>
<tbody>
<tr><th></th><th class="numcol">Бактэрыі</th><th class="numcol">Археі</th><th class="numcol">Усяго</th></tr>
<tr class="sumtab-row">
<td>Колькасць паслядоўнасцяў генаў</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Колькасць унікальных паслядоўнасцяў</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Колькасць відаў</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Колькасць геномаў</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Колькасць геномаў 1-й катэгорыі</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Колькасць геномаў 2-й катэгорыі</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Колькасць геномаў 3-й катэгорыі</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Total'] }}</td>
</tr>
</tbody>
</table>
<div id="gene-lens" class="pad-anchor"></div>
<table class="sum-table">
<caption>Даўжыня генаў 16S рРНК</caption>
<tbody>
<tr><th></th><th class="numcol">Бактэрыі</th><th class="numcol">Археі</th></tr>
<tr class="sumtab-row">
<td>Мінімальная даўжыня, п.а.</td>
<td class="numcol">{{ ribogrove_len_dict['min']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['min']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>25-ы працэнтыль, п.а. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['25perc']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['25perc']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Медыяна, п.а. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['median']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['median']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>75-ы працэнтыль, п.а. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['75perc']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['75perc']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Сярэдняя даўжыня, п.а. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['mean']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['mean']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Мода, п.а. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['modes']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['modes']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Максімальная даўжыня, п.а.</td>
<td class="numcol">{{ ribogrove_len_dict['max']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['max']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Сярэднеквадратычнае адхіленне, п.а.<sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['std']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['std']['Archaea'] }}</td>
</tr>
</tbody>
</table>
<p><sup>*</sup> Метрыкі, пазначаныя зорачкай, вылічаныя з папярэдняй нармалізацыяй: вылічэннем <em>медыян</em> даўжынь генаў у межах віду.</p>
<div id="copy-number" class="pad-anchor"></div>
<table class="sum-table">
<caption>Колькасць копій генаў 16S рРНК<br><small>(Бактэрыі і Археі)</small></caption>
<tbody>
<tr><th class="numcol">Колькасць копій <sup>*</sup></th><th class="numcol">Колькасць відаў</th><th class="numcol">Адсотак відаў, %</th></tr>
{% for _, row in ribogrove_copy_number_df.iterrows() %}
<tr class="sumtab-row">
<td class="numcol">{{ row['copy_number'] }}</td>
<td class="numcol">{{ row['number_of_species'] }}</td>
<td class="numcol">{{ row['percent_of_species'] }}</td>
</tr>
{% endfor %}
</tbody>
</table>
<p><sup>*</sup> Медыянная колькасць копій у межах віду.</p>

<div id="top-longest" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 найдаўжэйшых генаў 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Арганізм</th><th class="numcol">Даўжыня, п.а.</th><th class="alnleft">SeqID <sup>*</sup></th><th class="alnleft">Ідэнтыфікатар<br>геномнай зборкі</th></tr>
<tr>
<td colspan="4" class="subhead">Бактэрыі</td>
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
<td colspan="4" class="subhead">Археі</td>
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
<p><sup>*</sup> SeqID &mdash; гэта ідэнтыфікатар паслядоўнасці RiboGrove.</p>

<div id="top-shortest" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 найкарацейшых генаў 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Арганізм</th><th class="numcol">Даўжыня, п.а.</th><th class="alnleft">SeqID <sup>*</sup></th><th class="alnleft">Ідэнтыфікатар<br>геномнай зборкі</th></tr>
<tr>
<td colspan="4" class="subhead">Бактэрыі</td>
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
<td colspan="4" class="subhead">Археі</td>
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
<p><sup>*</sup> SeqID &mdash; гэта ідэнтыфікатар паслядоўнасці RiboGrove.</p>

<div id="top-copy-num" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 геномаў з найбольшай колькасцю копій генаў 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Арганізм</th><th class="numcol">Колькасць копій</th><th class="alnleft">Ідэнтыфікатар<br>геномнай зборкі</th></tr>
<tr>
<td colspan="4" class="subhead">Бактэрыі</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['title']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Археі</td>
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
<caption>Топ-10 геномаў з найбольшай зменлівасцю генаў 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Арганізм</th><th class="numcol">Сума энтрапіі<sup> *</sup>, біты</th><th class="numcol">Сярэдняя энтрапія<sup> *</sup>, біты</th><th class="numcol">Колькасць зменлівых пазіцый</th><th class="numcol">Колькасць копій генаў</th><th class="alnleft">Ідэнтыфікатар<br>геномнай зборкі</th></tr>
<tr>
<td colspan="6" class="subhead">Бактэрыі</td>
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
<td colspan="6" class="subhead">Археі</td>
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

<p><sup>*</sup> Энтрапія ў гэтым выпадку &mdash; гэта <a href="https://ru.wikipedia.org/wiki/Информационная_энтропия">энтрапія Шэнана</a>, якую вылічалі для кожнай калонкі множнага выраўноўвання поўнапамерных паслядоўнасцяў генаў 16S рРНК кожнага генома. Пасля гэтага падсумоўвалі энтрапію (калонка <q>Сума энтрапіі</q>) і усярэднівалі (калонка <q>Сярэдняя энтрапія</q>).</p>
<p><sup>**</sup> <i>Halomicrobium</i> sp. ZPS1 &mdash; паказальны выпадак. У гэтым геноме два гены 16S рРНК, таму энтрапія роўная колькасці нуклеатыдаў, якія не супадаюць у іхніх паслядоўнасцях. Вылічыўшы адсотак ідэнтычнасці (скарочана АІ) паслядоўнасцяў гэтых двух генаў, атрымаем 90,70%! Нагадаем, што за парог АІ паслядоўнасцяў генаў 16S рРНК для аднясення арганізмаў да адрозных <em>родаў</em> звычайна бярэцца (зрэшты, адвольна) значэнне 95%.</p>

<hr>
<div id="searching-data" class="pad-anchor"></div>
<h2>Пошук дадзеных у RiboGrove</h2>
<p>RiboGrove &mdash; вельмі мінімалістычная база дадзеных: яна складаецца з двух файлаў fasta з дадатковымі метададзенымі. Таму пашыраных сродкаў пошуку няма. Мы прызнае&#769;м гэтую праблему і пакідаем ніжэй некалькі парад наконт пошуку дадзеных у RiboGrove.</p>
<div id="header-format" class="pad-anchor"></div>
<h3>Фармат загалоўка</h3>
<p>Запісы паслядоўнасцяў RiboGrove у фармаце fasta маюць наступны фармат загалоўка:</p>
<p class="samp-highl samp-vwide">&gt;NZ_CP079719.1:86193-87742_plus Bacillus_velezensis ;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus; category:2</p>
<p>Асноўныя часткі загалоўка адасобленыя адна ад адной прабеламі. У загалоўку павінна быць чатыры такія часткі:</p>
<ol>
<li>Ідэнтыфікатар паслядоўнасці (SeqID): <span class="samp-highl">NZ_CP079719.1:86193-87742_plus</span>. SeqID, у сваю чаргу, складаецца з трох частак:<ol>
<li>код доступу паслядоўнасці RefSeq, з якой паходзіць паслядоўнасць гена: <span class="samp-highl">NZ_CP079719.1</span>;</li>
<li>каардынаты гена у гэтай геномнай паслядоўнасці: <span class="samp-highl">86193-87742</span> (нумерацыя ідзе з 1-кі, інтэрвалы уключаюць левыя і правыя межы);</li>
<li>ланцуг геномнай паслядоўнасці, на якой ген знаходзіцца: <span class="samp-highl">plus</span> (альбо <span class="samp-highl">minus</span>).</li>
</ol></li>
<li>Поўнае імя арганізма, згодна з <a href="https://ncbi.nlm.nih.gov/taxonomy">базай дадзеных NCBI Taxonomy</a>: <span class="samp-highl">Bacillus_velezensis</span>.</li>
<li>Радок таксаноміі, а менавіта назвы дамену (<span class="samp-highl">Bacteria</span>), тыпу (<span class="samp-highl">Firmicutes</span>), класу (<span class="samp-highl">Bacilli</span>), парадку (<span class="samp-highl">Bacillales</span>), сямейства (<span class="samp-highl">Bacillaceae</span>) ды роду (<span class="samp-highl">Bacillus</span>). Назвы адасобленыя і аточаныя кропкамі з коскай (<span class="samp-highl">;</span>).</li>
<li>Катэгорыя генома, з якога паслядоўнасць гена паходзіць: (<span class="samp-highl">category:2</span>). Геномы катэгорыі 1 найбольш надзейныя, а катэгорыі 3 &mdash; найменш надзейныя.</li>
</ol>
<h3>Адбор паслядоўнасцяў</h3>
<p>Адбіраць тыя ці іншыя паслядоўнасці з файлаў fasta можна з дапамогай праграмы <strong>Seqkit</strong> (<a href="https://github.com/shenwei356/seqkit">ГітХаб-рэпазіторый</a>, <a href="https://bioinf.shenwei.me/seqkit/">дакументацыя</a>). Гэта бясплатная, кросплатформавая, шматфункцыянальная і даволі шустрая праграма, якая можа апрацоўваць як &bdquo;гзіп&rsquo;нутыя&ldquo;, так і простыя нясціснутыя файлы fasta. Для адбору паслядоўнасцяў зручна выкарыстоўваць падпраграмы <span class="samp">seqkit grep</span> і <span class="samp">seqkit seq</span> .</p>
<div id="subset-by-header" class="pad-anchor"></div>
<h4>Адбор паслядоўнасцяў паводле загалоўкаў</h4>
<p>Няхай мы робім выбарку з файла <span class="samp">ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</span>. Карыснымі могуць быць наступныя прыклады каманд карыстання праграмай <span class="samp">seqkit grep</span>:</p>
<p><strong>Прыклад 1</strong>. Выбраць адну паслядоўнасць паводле SeqID.</p>
<p class="samp-highl samp-vwide">seqkit grep -p "NZ_CP079719.1:86193-87742_plus" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опцыя <span class="samp">-p</span> задае радок, які праграма будзе шукаць у fasta-загалоўках (на самой справе, толькі ў той частцы загалоўкаў, дзе запісаны SeqID: да першага прабела).</p>
<p><strong>Прыклад 2</strong>. Выбраць усе паслядоўнасці генаў адной геномнай паслядоўнасці RefSeq паводле кода доступу.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "NZ_CP079719.1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Тут спатрэбяцца яшчэ дзве опцыі: <span class="samp">-n</span> і <span class="samp">-r</span>. Першая кажа праграме шукаць падрадкі&#769; у загалоўках цалкам, а не толькі ў ідэнтыфікатарах (SeqID). Другая опцыя кажа праграме, што ёй трэба шукаць не толькі тыя загалоўкі, якія цалкам супадаюць з шуканым радком, а таксама і тыя, якія змяшчаюць шуканы радок як сваю частку.</p>
<p><strong>Прыклад 3</strong>. Выбраць усе паслядоўнасці генаў аднаго генома (ідэнтыфікатар геномнай зборкі <a href="https://ncbi.nlm.nih.gov/assembly/10577151">10577151</a>), які мае два рэпліконы: <a href="https://ncbi.nlm.nih.gov/nuccore/NZ_CP079110.1">NZ_CP079110.1</a> і <a href="https://ncbi.nlm.nih.gov/nuccore/NZ_CP079111.1">NZ_CP079111.1</a>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nr -p "NZ_CP079110.1" -p "NZ_CP079111.1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Прыклад 4</strong>. Выбраць усе паслядоўнасці актынабактэрый.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ";Actinobacteria;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">На ўсялякі выпадак, атачайце таксаномію кропкамі з коскай (<span class="samp-highl">;</span>).</p>
<p><strong>Прыклад 5</strong>. Выбраць усе паслядоўнасці з геномаў катэгорыі 1.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "category:1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Прыклад 6</strong>. Выбраць усе паслядоўнасці апроч тых, якія належаць да тыпу <i>Firmicutes</i>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nvrp ";Firmicutes;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Звярніце ўвагу на опцыю <span class="samp">-v</span> у шэгару опцый <span class="samp">-nvrp</span>. Гэтая опцыя інвертуе паводзіны праграмы, і цяпер яна знойдзе ўсе паслядоўнасці, загалоўкі якіх <em>не</em> змяшчаюць радок <span class="samp">;Firmicutes;</span>.</p>
<div id="subset-by-len" class="pad-anchor"></div>
<h4>Адбор паслядоўнасцяў паводле даўжыні</h4>
<p>Адбіраць паслядоўнасці паводле даўжыні можна з дапамогай праграмы <span class="samp">seqkit seq</span>.</p>
<p><strong>Прыклад 1</strong>. Выбраць усе паслядоўнасці даўжэй за 1600 п.а.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1601 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опцыя <span class="samp">-m</span> задае мінімальную даўжыню паслядоўнасцяў, якія праграма падасць на выхад.</p>
<p><strong>Прыклад 2</strong>. Выбраць усе паслядоўнасці карацей за 1500 п.а.</p>
<p class="samp-highl samp-vwide">seqkit seq -M 1499 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опцыя <span class="samp">-M</span> задае максімальную даўжыню паслядоўнасці, якія праграма падасць на выхад.</p>
<p><strong>Прыклад 3</strong>. Выбраць усе паслядоўнасці, даўжыня якіх знаходзіцца ў прамежку [1500, 1600] п.а.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1500 -M 1600 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<div id="select-head" class="pad-anchor"></div>
<h4>Адбор дадзеных загалоўкаў</h4>
<p>Часам бывае карысна атрымаць толькі дадзеныя загалоўкаў з файла fasta. Гэта можна зрабіць з дапамогай праграмы <span class="samp">seqkit seq</span>.</p>
<p><strong>Прыклад 1</strong>. Выбраць проста ўсе загалоўкі.</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опцыя <span class="samp">-n</span> кажа праграме, каб тая падавала на выхад толькі загалоўкі.</p>
<p><strong>Прыклад 2</strong>. Выбраць усе ідэнтыфікатары паслядоўнасцяў (то бок часткі загалоўкаў да першага прабела).</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опцыя <span class="samp">-i</span> кажа праграме, каб тая падавала на выхад толькі ідэнтыфікатары.</p>
<p><strong>Прыклад 3</strong>. Выбраць усе (RefSeq&rsquo;аўскія) коды доступу.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f1 -d':' | sort | uniq</p>
<p class="samp-comment">Такая команда спрацуе, толькі калі на камп&rsquo;ютары ўсталяваныя ўтыліты <span class="samp">sort</span>, <span class="samp">cut</span> і <span class="samp">uniq</span> (на Linux і Mac OS гэтыя праграмы звычайна ўсталяваныя разам з аперацыйная сістэмай).</p>
<br>
<hr>
<p>RiboGrove, {{ ribogrove_release_date }}</p>
</div>
</body>
</html>
