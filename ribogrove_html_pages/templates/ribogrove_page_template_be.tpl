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
<p>&#127760;&nbsp;<a href="index.html">English</a> / <span class="active-lang">Беларуская</span> / <a href="latest_ribogrove_release_ua.html">Українська</a> / <a href="latest_ribogrove_release_ru.html">Русский</a></p>
</div>
<p><div class="warn-section"><a href="https://mbio.bas-net.by/cager/en/ribogrove">Асноўны сайт</a>, дзе размяшчаецца база дадзеных RiboGrove, можа быць недасяжным па-за межамі Беларусі праз тэхнічныя прычыны ды бяду даогул. Таму існуе гэтае люстэрка, а файлы RiboGrove дасяжныя праз Дропбокс (спасылкі <a href="#downloads">ніжэй</a>).</div></p><br>
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
<li><a href="#downloads">Файлы</a>
<ul class="ribogrove-nested-list">
<li><a href="#current-release">Актуальны выпуск RiboGrove — {{ ribogrove_release_number }}</li>
<li><a href="#release-archive">Архіў выпускаў RiboGrove</li>
<li><a href="#release-notes">Каментары да выпуску</li>
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
<li><a href="#primers-coverages">Спектр дзеяння пар праймераў да розных V-регіёнаў бактэрыяльных генаў 16S рРНК</a></li>
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
<li><a href="#contacts">Кантакты</a></li>
<li><a href="#citing-ribogrove">Цытаванне RiboGrove</a></li>
<li><a href="#faq">Пытанні, якія ў людзей узнікаюць пра RiboGrove</a></li>
</ul>
<hr>
<div id="overview" class="pad-anchor"></div>
<h2>Што такое RiboGrove</h2>
<p>RiboGrove — гэта база дадзеных паслядоўнасцей генаў 16S рРНК бактэрый і архей.</p>
<p>RiboGrove асноўваецца на базе дадзеных <a href="https://ncbi.nlm.nih.gov/refseq/">RefSeq</a>. RiboGrove змяшчае толькі поўнапамерныя паслядоўнасці генаў 16S рРНК, а самі паслядоўнасці паходзяць з цалкам сабраных геномаў прокарыёт, дэпанаваных у RefSeq. Таму мы пастулюем высокую надзейнасць паслядоўнасцей, размешчаных у RiboGrove.</p>
<div id="other-databases" class="pad-anchor"></div>
<p>Табліца ніжэй паказвае, чым (якасна) RiboGrove адрозніваецца ад падобных баз дадзеных паслядоўнасцей РНК, а менавіта <a href="https://rrndb.umms.med.umich.edu/"><i>rrn</i>DB</a>, <a href="https://www.arb-silva.de/">Silva</a>, <a href="http://rdp.cme.msu.edu/">RDP</a> і <a href="https://greengenes.secondgenome.com/">Greengenes</a>. Кажучы каротка, RiboGrove саступае аналагам у колькасці паслядоўнасцей і іх разнастайнасці, але пераўзыходзіць у надзейнасці паслядоўнасцей.</p>
<table class="sum-table">
<tbody>
  <tr class="sumtab-row"><th></th><th class="alcentr">RiboGrove</th><th class="alcentr"><i>rrn</i>DB</th><th class="alcentr">Silva</th><th class="alcentr">RDP</th><th class="alcentr">Greengenes</th></tr>
  <tr class="sumtab-row"><td class="feature-col">Прадстаўленыя арганізмы</td><td class="alcentr">Бактэрыі<br>Археі</td><td class="alcentr">Бактэрыі<br>Археі</td><td class="alcentr">Бактэрыі<br>Археі<br>Эўкарыёты</td><td class="alcentr">Бактэрыі<br>Археі<br>Эўкарыёты</td><td class="alcentr">Бактэрыі<br>Археі</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Прадстаўленыя субадзінкі рыбасом</td><td class="alcentr">Малая</td><td class="alcentr">Малая</td><td class="alcentr">Вялікая<br>Малая</td><td class="alcentr">Вялікая<br>Малая</td><td class="alcentr">Малая</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Змяшчае паслядоўнасці сабраных геномаў</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Змяшчае паслядоўнасці ампліконаў</td><td class="alcentr">Не</td><td class="alcentr">Не</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Змяшчае няпоўныя паслядоўнасці генаў</td><td class="alcentr">Не</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Адрознівае катэгорыі геномаў</td><td class="alcentr">Так</td><td class="alcentr">Не</td><td class="alcentr"><i>Крытэр<br>непрымянімы</i></td><td class="alcentr"><i>Крытэр<br>непрымянімы</i></td><td class="alcentr"><i>Крытэр<br>непрымянімы</i></td></tr>
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
  <li>Зборка геному утрымлівае хаця б адзін запіс базы дадзеных RefSeq, у назве якой ёсць фраза «map unlocalized», і гэты запіс утрымлівае паслядоўнасць гена 16S рРНК або яго частку.</li>
</ul>
<br>
<p>Праграмы, з дапамогай якіх была створана RiboGrove, знаходзяцца ў наступным ГітХаб-рэпазіторыі: <a href="https://github.com/masikol/ribogrove-tools">ribogrove-tools</a>.</p>
<hr>
<div id="downloads" class="pad-anchor"></div>
<h2>Файлы</h2>
<div id="current-release" class="pad-anchor"></div>
<h3>Актуальны выпуск RiboGrove — {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
<p>Выпуск заснованы на базе дадзеных RefSeq {{ refseq_release }}.</p>
<ul>
<li>Fasta-файл з поўнапамернымі паслядоўнасцямі генаў 16S рРНК. <a href="ribogrove_releases/{{ ribogrove_release_number }}/ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Спампаваць („гзіпʼнуты“ fasta-файл, {{ final_fasta_fsize_fmt }} Мб)</a></li>
<li>Метададзеныя. <a href="ribogrove_releases/{{ ribogrove_release_number }}/metadata_ribogrove_{{ ribogrove_release_number }}.zip">Спампаваць (zip-архіў {{ metadata_fsize_fmt }} Мб)</a><br>
<details>
  <summary><small>Якая менавіта інфармацыя утрымліваецца ў метададзеных?</small></summary>
  <small>
    <p>Метададзеныя складаюцца з наступных файлаў:</p>
    <ol>
    <li>
      <span class="samp-highl">discarded_sequences.fasta.gz</span><br>
      Гэта fasta-файл, які змяшчае паслядоўнасці, якія прысутнічалі ў зыходных геномах, былі пры гэтым анатаваныя як гены 16S рРНК, але якія былі „забракаваныя“ (discarded) праз іхнюю непаўнату, прысутнасць унутрыгенных паўтораў і г.д. Такім чынам, гэтыя паслядоўнасці не ўключаныя ў RiboGrove.
    </li>
    <li>
      <span class="samp-highl">source_RefSeq_genomes.tsv</span><br>
      Гэта TSV-файл, які змяшчае інфармацыю пра тое, якія геномы выкарыстоўваліся для стварэння RiboGrove.
    </li>
    <li>
      <span class="samp-highl">gene_seqs_statistics.tsv</span>, <span class="samp-highl">discarded_gene_seqs_statistics.tsv</span><br>
      Гэта TSV-файлы, якія змяшчаюць інфармацыю пра нуклеатыдны склад, памер, геномную і таксанамічную прыналежнасць паслядоўнасцей генаў. Першы файл з гэтых двух апісвае фінальныя паслядоўнасці, уключаныя ў RiboGrove, а другі — „забракаваныя“ паслядоўнасці.
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
<p>Файл fasta сціснуты з дапамогаю праграмы gzip, а файлы з метададзенымі сціснутыя ў zip-архіў. Каб распакаваць іх, карыстальнікі Linux і Mac OS могуць скарастыцца якраз праграмамі gzip і zip; гэтыя праграмы звычайна ўсталяваныя разам з аперыцыйнай сістэмай. Карыстальнікам Windows можа прыдацца бясплатная праграма-архіватар <a href="https://www.7-zip.org/">7-Zip</a>.</p>
<div id="release-archive" class="pad-anchor"></div>
<h3>Архіў выпускаў RiboGrove</h3>
<p>Усе папярэднія выпускі можна знайсці ў <a href="archive/ribogrove_archive_en.html"> архіве выпускаў RiboGrove</a>. На жаль, ён вядзецца толькі па-англійску.</p>
<div id="release-notes" class="pad-anchor"></div>
<h3>Каментары да выпуску</h3>
<p><i>Няма прынцыповых адрозненняў ад папярэдняга выпуску.</i></p>
<p>Каментары да ўсіх выпускаў RiboGrove можна знайсці на <a href="ribogrove_release_notes_en.html">старонцы каментароў да выпускаў</a>.</p>
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
<caption>Колькасць копій генаў 16S рРНК</caption>
<tbody>
<tr><th class="numcol">Колькасць копій <sup>*</sup></th><th class="alncenter" colspan="2">Бактэрыі</th><th class="alncenter" colspan="2">Археі</th></tr>
<tr><td></td><td>Колькасць відаў</td><td>Адсотак відаў, %</td><td>Колькасць відаў</td><td>Адсотак відаў, %</td></tr>
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
<p><sup>*</sup> Медыянная колькасць копій у межах віду.</p>

<div id="top-longest" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 найдаўжэйшых генаў 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Арганізм</th><th class="numcol">Даўжыня, п.а.</th><th class="alnleft">SeqID <sup>*</sup></th><th class="alnleft">Код доступу<br>геномнай зборкі</th></tr>
<tr>
<td colspan="4" class="subhead">Бактэрыі</td>
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
<td colspan="4" class="subhead">Археі</td>
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
<p><sup>*</sup> SeqID — гэта ідэнтыфікатар паслядоўнасці RiboGrove.</p>

<div id="top-shortest" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 найкарацейшых генаў 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Арганізм</th><th class="numcol">Даўжыня, п.а.</th><th class="alnleft">SeqID <sup>*</sup></th><th class="alnleft">Код доступу<br>геномнай зборкі</th></tr>
<tr>
<td colspan="4" class="subhead">Бактэрыі</td>
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
<td colspan="4" class="subhead">Археі</td>
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
<p><sup>*</sup> SeqID — гэта ідэнтыфікатар паслядоўнасці RiboGrove.</p>

<div id="top-copy-num" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 геномаў з найбольшай колькасцю копій генаў 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Арганізм</th><th class="numcol">Колькасць копій</th><th class="alnleft">Код доступу<br>геномнай зборкі</th></tr>
<tr>
<td colspan="4" class="subhead">Бактэрыі</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Археі</td>
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
<table class="sum-table">
<caption>Топ-10 геномаў з найбольшай зменлівасцю генаў 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Арганізм</th><th class="numcol">Сума энтрапіі<sup> *</sup>, біты</th><th class="numcol">Сярэдняя энтрапія<sup> *</sup>, біты</th><th class="numcol">Колькасць зменлівых пазіцый</th><th class="numcol">Колькасць копій генаў</th><th class="alnleft">Код доступу<br>геномнай зборкі</th></tr>
<tr>
<td colspan="6" class="subhead">Бактэрыі</td>
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
<td colspan="6" class="subhead">Археі</td>
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

<p><sup>*</sup> Энтрапія ў гэтым выпадку — гэта <a href="https://ru.wikipedia.org/wiki/Информационная_энтропия">энтрапія Шэнана</a>, якую вылічалі для кожнай калонкі множнага выраўноўвання поўнапамерных паслядоўнасцяў генаў 16S рРНК кожнага геному. Пасля гэтага падсумоўвалі энтрапію (калонка «Сума энтрапіі») і усярэднівалі (калонка «Сярэдняя энтрапія»).</p>
<p><sup>**</sup> <i>Halomicrobium</i> sp. ZPS1 — паказальны выпадак. У гэтым геноме два гены 16S рРНК, таму энтрапія роўная колькасці нуклеатыдаў, якія не супадаюць у іхніх паслядоўнасцях. Вылічыўшы адсотак ідэнтычнасці (скарочана АІ) паслядоўнасцяў гэтых двух генаў, атрымаем 90,70%! Нагадаем, што за парог АІ паслядоўнасцяў генаў 16S рРНК для аднясення арганізмаў да адрозных <em>родаў</em> звычайна бярэцца (зрэшты, адвольна) значэнне 95%.</p>

<div id="primers-coverages" class="pad-anchor"></div>
<table class="sum-table"><caption>Спектр дзеяння<sup>*</sup> пар праймераў да розных V-регіёнаў бактэрыяльных генаў 16S рРНК</caption>
<tbody class="primer-cov-tbody">
<tr>
  <th class="alnleft" rowspan="2">Тып<br />бактэрый</th>
  <th class="numcol" rowspan="2">Колькасць<br />геномаў</th>
  <th class="numcol">Поўны<br />ген</th>
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

<p><sup>*</sup> Спектр дзеяння пары праймераў — гэта адсотак геномаў, якія маюць прынамсі адзін ген 16S рРНК, які магчыма ампліфікаваць з дапамогай ПЦР і з выкарыстаннем гэтай пары праймераў. Падрабязнасці апісаныя ў нашым <a href="https://masikol.github.io/latest_ribogrove_release_be.html#citing-ribogrove">артыкуле пра RiboGrove</a>.</p>

<table class="sum-table"><caption>Праймеры, спектр дзеяння якіх ацэньваўся</caption>
<tbody>
  <tr><th>Назва праймера</th><th>Паслядоўнасць</th><th>Спасылка</th></tr>
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

<p><sup>*</sup> Праймеры 341F і 785R выкарыстоўваюцца ў <a href="https://support.illumina.com/downloads/16s_metagenomic_sequencing_library_preparation.html">пратаколе</a> падрыхтоўкі бібліятэкі для секвенавання ампліконаў регіёнаў V3–V4 генаў 16S рРНК на апараце Illumina MiSeq.</p>

<hr>
<div id="searching-data" class="pad-anchor"></div>
<h2>Пошук дадзеных у RiboGrove</h2>
<p>RiboGrove — вельмі мінімалістычная база дадзеных: яна складаецца з двух файлаў fasta з дадатковымі метададзенымі. Таму пашыраных сродкаў пошуку няма. Мы прызнае&#769;м гэтую праблему і пакідаем ніжэй некалькі парад наконт пошуку дадзеных у RiboGrove.</p>
<div id="header-format" class="pad-anchor"></div>
<h3>Фармат загалоўка</h3>
<p>Запісы паслядоўнасцяў RiboGrove у фармаце fasta маюць наступны фармат загалоўка:</p>
<p class="samp-highl samp-vwide" style="font-size:14px">&gt;GCF_000978375.1:NZ_CP009686.1:8908-10459:plus ;d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__cereus; category:1</p>
<p>Асноўныя часткі загалоўка адасобленыя адна ад адной прабеламі. У загалоўку павінна быць тры такія часткі:</p>
<ol>
<li>Ідэнтыфікатар паслядоўнасці (SeqID): <span class="samp-highl">GCF_000978375.1:NZ_CP009686.1:8908-10459:plus</span>. SeqID, у сваю чаргу, складаецца з чатырох частак, адасобленых двукропʼямі (<span class="samp-highl">:</span>):
<ol>
<li>Код доступу зборкі геному, з якога паходзіць паслядоўнасць гена: <span class="samp-highl">GCF_000978375.1</span>.</li>
<li>Код доступу паслядоўнасці RefSeq, з якой паходзіць паслядоўнасць гена: <span class="samp-highl">NZ_CP009686.1</span>.</li>
<li>Каардынаты гена у гэтай геномнай паслядоўнасці: <span class="samp-highl">8908-10459</span> (нумерацыя ідзе з 1-кі, інтэрвалы уключаюць левыя і правыя межы).</li>
<li>Ланцуг геномнай паслядоўнасці, на якой ген знаходзіцца: <span class="samp-highl">plus</span> (альбо <span class="samp-highl">minus</span>).</li>
</ol>
</li>
<li>Радок таксаноміі, а менавіта назвы дамену (<span class="samp-highl">Bacteria</span>), тыпу (<span class="samp-highl">Firmicutes</span>), класу (<span class="samp-highl">Bacilli</span>), парадку (<span class="samp-highl">Bacillales</span>), сямейства (<span class="samp-highl">Bacillaceae</span>), роду (<span class="samp-highl">Bacillus</span>), а таксама відавы эпітэт (<span class="samp-highl">cereus</span>).<br>
Перад кожнай назвай пішацца прэфікс, які пазначае ранг: <span class="samp-highl">d__</span> для дамену, <span class="samp-highl">p__</span> для тыпу, <span class="samp-highl">c__</span> для класу, <span class="samp-highl">o__</span> для парадку, <span class="samp-highl">f__</span> для сямейства, <span class="samp-highl">g__</span> для роду і <span class="samp-highl">s__</span> для відавага эпітэту. У прэфіксах пасля літары пішуцца <strong>два</strong> знакі ніжняга падкрэслівання.<br>
Таксанамічныя назвы адасобленыя і аточаныя кропкамі з коскай (<span class="samp-highl">;</span>).</li>
<li>Катэгорыя геному, з якога паслядоўнасць гена паходзіць: (<span class="samp-highl">category:1</span>).</li>
</ol>
<h3>Адбор паслядоўнасцяў</h3>
<p>Адбіраць тыя ці іншыя паслядоўнасці з файлаў fasta можна з дапамогай праграмы <strong>Seqkit</strong> (<a href="https://github.com/shenwei356/seqkit">ГітХаб-рэпазіторый</a>, <a href="https://bioinf.shenwei.me/seqkit/">дакументацыя</a>). Гэта бясплатная, кросплатформавая, шматфункцыянальная і даволі шустрая праграма, якая можа апрацоўваць як „гзіпʼнутыя“, так і простыя нясціснутыя файлы fasta. Для адбору паслядоўнасцяў зручна выкарыстоўваць падпраграмы <span class="samp">seqkit grep</span> і <span class="samp">seqkit seq</span> .</p>
<div id="subset-by-header" class="pad-anchor"></div>
<h4>Адбор паслядоўнасцяў паводле загалоўкаў</h4>
<p>Няхай мы робім выбарку з файла <span class="samp">ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</span>. Карыснымі могуць быць наступныя прыклады каманд карыстання праграмай <span class="samp">seqkit grep</span>:</p>
<p><strong>Прыклад 1</strong>. Выбраць адну паслядоўнасць паводле SeqID.</p>
<p class="samp-highl samp-vwide">seqkit grep -p "GCF_000978375.1:NZ_CP009686.1:8908-10459:plus" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опцыя <span class="samp">-p</span> задае радок, які праграма будзе шукаць у fasta-загалоўках (на самой справе, толькі ў той частцы загалоўкаў, дзе запісаны SeqID: да першага прабела).</p>
<p><strong>Прыклад 2</strong>. Выбраць усе паслядоўнасці генаў адной геномнай паслядоўнасці RefSeq паводле кода доступу <span class="samp">NZ_CP009686.1</span>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ":NZ_CP009686.1:" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Тут спатрэбяцца яшчэ дзве опцыі: <span class="samp">-n</span> і <span class="samp">-r</span>. Першая кажа праграме шукаць падрадкі&#769; у загалоўках цалкам, а не толькі ў ідэнтыфікатарах (SeqID). Другая опцыя кажа праграме, што ёй трэба шукаць не толькі тыя загалоўкі, якія цалкам супадаюць з шуканым радком, а таксама і тыя, якія змяшчаюць шуканы радок як сваю частку.</p>
<p class="samp-comment">Каб забяспечыць спецыфічнасць пошуку, атачайце код доступу двукропʼямі (<span class="samp-highl">:</span>).</p>
<p><strong>Прыклад 3</strong>. Выбраць усе паслядоўнасці генаў аднаго геному (код доступу геномнай зборкі <a href="https://ncbi.nlm.nih.gov/assembly/GCF_019357495.1">GCF_019357495.1</a>).</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "GCF_019357495.1:" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Каб забяспечыць спецыфічнасць пошуку, пішыце двукропʼе (<span class="samp-highl">:</span>) пасля коду доступу зборкі.</p>
<p><strong>Прыклад 4</strong>. Выбраць усе паслядоўнасці актынабактэрый.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ";p__Actinobacteria;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Каб забяспечыць спецыфічнасць пошуку, атачайце таксаномію кропкамі з коскай (<span class="samp-highl">;</span>).</p>
<p><strong>Прыклад 5</strong>. Выбраць усе паслядоўнасці з геномаў катэгорыі 1.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "category:1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Прыклад 6</strong>. Выбраць усе паслядоўнасці апроч тых, якія належаць да тыпу <i>Firmicutes</i>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nvrp ";p__Firmicutes;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Звярніце ўвагу на опцыю <span class="samp">-v</span> у шэгару опцый <span class="samp">-nvrp</span>. Гэтая опцыя інвертуе паводзіны праграмы, і цяпер яна знойдзе ўсе паслядоўнасці, у загалоўках якіх <em>няма</em> радку «<span class="samp">;p__Firmicutes;</span>».</p>
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
<p><strong>Прыклад 3</strong>. Выбраць усе коды доступу геномных паслядоўнасцяў з БД RefSeq.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f2 -d':' | sort | uniq</p>
<p class="samp-comment">Такая команда спрацуе, толькі калі на кампʼютары ўсталяваныя ўтыліты <span class="samp">cut</span>, <span class="samp">sort</span>, і <span class="samp">uniq</span> (на Linux і Mac OS гэтыя праграмы звычайна ўсталяваныя разам з аперацыйная сістэмай).</p>
<p><strong>Прыклад 4</strong>. Выбраць усе коды доступу геномных зборак.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f1 -d':' | sort | uniq</p>
<p class="samp-comment">Такая команда спрацуе, толькі калі на кампʼютары ўсталяваныя ўтыліты <span class="samp">cut</span>, <span class="samp">sort</span> і <span class="samp">uniq</span> (на Linux і Mac OS гэтыя праграмы звычайна ўсталяваныя разам з аперацыйная сістэмай).</p>
<p><strong>Прыклад 5</strong>. Выбраць усе назвы тыпаў арганізмаў.</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | grep -Eo ';p__[^;]+' | sed -E 's/;|p__//g' | sort | uniq</p>
<p class="samp-comment">Такая команда спрацуе, толькі калі на кампʼютары ўсталяваныя ўтыліты <span class="samp">grep</span>, <span class="samp">sed</span>, <span class="samp">sort</span> і <span class="samp">uniq</span> (на Linux і Mac OS гэтыя праграмы звычайна ўсталяваныя разам з аперацыйная сістэмай).</p>
<hr>
<div id="contacts" class="pad-anchor"></div>
<h2>Кантакты</h2>
<p>Калі вы маеце пытанне пра RiboGrove, калі ласка, звяртайцеся да Максіма Сікаленкі на эл. адрас sikolenko<img class="sabaka" src="img/pes2.gif" alt="[ at ]" align="bottom">bio.bsu.by.</p>
<hr>
<div id="citing-ribogrove" class="pad-anchor"></div>
<h2>Цытаванне RiboGrove</h2>
<p>Калі база дадзеных RiboGrove прыдалася вам у вашых даследаваннях, спасылайцеся, калі ласка, на наступную публікацыю:</p>
<p>Maxim A. Sikolenko, Leonid N. Valentovich. “RiboGrove: a database of full-length prokaryotic 16S rRNA genes derived from completely assembled genomes” // Research in Microbiology, Volume 173, Issue 4, May 2022, 103936.<br> (<span style="cursor:help"title="Лічбавы ідэнтыфікатар абʼекта">ЛІА</span>: <a href="https://doi.org/10.1016/j.resmic.2022.103936">10.1016/j.resmic.2022.103936</a>).</p>
<hr>
<div id="faq" class="pad-anchor"></div>
<h2>Пытанні, якія ў людзей узнікаюць пра RiboGrove</h2>
<div><details>
  <summary><b>1. Як зрабіць файл з таксаноміяй, з якім маглі бы працаваць праграмы з пакета QIIME2?</b></summary>
  <p>Для гэтага можаце скарыстацца скрыптом <a href="https://github.com/masikol/ribogrove-tools/blob/main/data_transformation_scripts/make_qiime_taxonomy_file.py"><span class="samp">make_qiime_taxonomy_file.py</span></a>: ён трансфармуе файл <span class="samp">metadata/taxonomy.tsv</span> у файл, з якім змогуць працаваць праграмы QIIME2. Даведку-кіраўніцтва для гэтага скрыпта можаце знайсці ў адпаведным <a href="https://github.com/masikol/ribogrove-tools/blob/main/data_transformation_scripts/README.md">файле README</a>.</p>
</details></div><br>
<div><details>
  <summary><b>2. Як захаваць у файл дадзеныя, адабраныя з дапамогай Seqkit?</b></summary>
  <p>У адпаведнай дыскусіі на форуме ўжо далі некалькі карысных адказаў і парад: <a href="https://bioinformatics.stackexchange.com/questions/20915/how-do-i-save-selected-sequences-in-seqkit-to-a-file">https://bioinformatics.stackexchange.com/questions/20915/how-do-i-save-selected-sequences-in-seqkit-to-a-file</a>.</p>
</details></div><br>
<div><details>
  <summary><b>3. Як знайсці патрэбныя запісы у fasta-файле паводле паслядоўнасці з дапамогай Seqkit?</b></summary>
  <p>У адпаведнай дыскусіі на форуме ўжо далі некалькі карысных адказаў і парад: <a href="https://www.biostars.org/p/9561418">https://www.biostars.org/p/9561418</a>.</p>
</details></div>
<br>
<hr>
<p>RiboGrove, {{ ribogrove_release_date }}</p>
</div>
</body>
</html>
