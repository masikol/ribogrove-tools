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
<p>&#127760;&nbsp;<a href="index.html">English</a> / <a href="latest_ribogrove_release_be.html">Беларуская</a> / <span class="active-lang">Українська</span> / <a href="latest_ribogrove_release_ru.html">Русский</a></p>
</div>
<p><div class="warn-section"><a href="https://mbio.bas-net.by/cager/en/ribogrove">Головний сайт</a>, де ми розміщуємо базу даних RiboGrove, може бути недосяжним поза межами Білорусі через технічні проблеми та біду взагалі. Актуальний реліз RiboGrove можна знайти тут — нехай база даних буде досяжна принаймни такою мірою.</div></p><br>
<h1 style="display:none">RiboGrove</h1>
<img src="img/RiboGrove_logo.png"><br><br>
<h2>Зміст</h2>
<ul>
{% if not archive %}
<li><a href="#overview">Що таке RiboGrove</a>
<ul class="ribogrove-nested-list">
<li><a href="#other-databases">RiboGrove і інші бази даних 16S рРНК</a></li>
<li><a href="#genome-categories">Категорії геномів</a></li>
</ul>
</li>
{% endif %}
<li><a href="#downloads">Файли</a>
<ul class="ribogrove-nested-list">
  {% if not archive %}
  <li><a href="#current-release">Актуальний реліз RiboGrove — {{ ribogrove_release_number }}</li>
  {% else %}
  <li><a href="#current-release">Реліз RiboGrove {{ ribogrove_release_number }}</li>
  {% endif %}
  <li><a href="#release-archive">Архів релізів RiboGrove</li>
  <li><a href="#release-notes">Коментарі до релізу</li>
</ul>
</li>
<li><a href="#summary">Статистичний опис</a>
<ul class="ribogrove-nested-list">
<li><a href="#db-size">Розмір RiboGrove</a></li>
<li><a href="#gene-lens">Довжина генів 16S рРНК</a></li>
<li><a href="#copy-number">Кількість копій генів 16S рРНК</a></li>
<li><a href="#top-longest">Топ-10 найдовших генів 16S рРНК</a></li>
<li><a href="#top-shortest">Топ-10 найкоротших генів 16S рРНК</a></li>
<li><a href="#top-copy-num">Топ-10 геномів за найбільшею кілкістю копій генів 16S рРНК</a></li>
<li><a href="#top-var">Топ-10 геномів за найбільшею мінливістю генів 16S рРНК</a></li>
<li><a href="#primers-coverages">Спектр диї пар праймерів до різних V-регіонів бактеріальних генів 16S рРНК</a></li>
</ul>
</li>
<li><a href="#searching-data">Пошук даних у RiboGrove</a>
<ul class="ribogrove-nested-list">
  <li><a href="#header-format">Формат заголовку</a></li>
<li><a href="#subset-by-header">Вибір послідовностей по заголовку</a></li>
<li><a href="#subset-by-len">Вибір послідовностей по довжині</a></li>
<li><a href="#select-head">Вибір даних заголовків</a></li>
</ul>
</li>
{% if not archive %}
<li><a href="#contacts">Контакти</a></li>
<li><a href="#citing-ribogrove">Цитування RiboGrove</a></li>
<li><a href="#faq">Питання, які люди задають про RiboGrove</a></li>
{% endif %}
</ul>
<hr>
{% if not archive %}
<div id="overview" class="pad-anchor"></div>
<h2>Що таке RiboGrove</h2>
<p>RiboGrove — це база даних послідовностей генів 16S рРНК бактерій та архей.</p>
<p>RiboGrove базується на базі даних <a href="https://ncbi.nlm.nih.gov/refseq/">RefSeq</a>. RiboGrove утримує тільки повнорозмірні послідовності генів 16S рРНК, а самі послідовності походять із повністю зібраних геномів прокаріотів, депонуваних у RefSeq. Тому ми постулюємо високу надійність послідовностей, розміщених у RiboGrove.</p>
<div id="other-databases" class="pad-anchor"></div>
<h3>RiboGrove і інші бази даних 16S рРНК</h3>
<p>Таблиця нижче показує, чим (якісно) RiboGrove відрізнювається від подібних баз даних послідовностей РНК, а саме <a href="https://rrndb.umms.med.umich.edu/"><i>rrn</i>DB</a>, <a href="https://www.arb-silva.de/">Silva</a>, <a href="http://rdp.cme.msu.edu/">RDP</a> та <a href="https://greengenes.secondgenome.com/">Greengenes</a>. Кажучи коротко, RiboGrove менша за аналоги щодо кількості послідовностей та їх різноманітності, але виграє у надійності послідовностей.</p>
<table class="sum-table">
<tbody>
  <tr class="sumtab-row"><th></th><th class="alcentr">RiboGrove</th><th class="alcentr"><i>rrn</i>DB</th><th class="alcentr">Silva</th><th class="alcentr">RDP</th><th class="alcentr">Greengenes</th></tr>
  <tr class="sumtab-row"><td class="feature-col">Репрезентовані організми</td><td class="alcentr">Бактерії<br>Археї</td><td class="alcentr">Бактерії<br>Археї</td><td class="alcentr">Бактерії<br>Археї<br>Еукаріоти</td><td class="alcentr">Бактерії<br>Археї<br>Еукаріоти</td><td class="alcentr">Бактерії<br>Археї</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Репрезентовані субодиниці рибосом</td><td class="alcentr">Мала</td><td class="alcentr">Мала</td><td class="alcentr">Велика<br>Мала</td><td class="alcentr">Велика<br>Мала</td><td class="alcentr">Мала</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Утримує послідовності зібраних геномів</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Утримує послідовності ампліконів</td><td class="alcentr">Ні</td><td class="alcentr">Ні</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Утримує неповні послідовності генів</td><td class="alcentr">Ні</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td><td class="alcentr">Так</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Відрізняє категориї геномів</td><td class="alcentr">Так</td><td class="alcentr">Ні</td><td class="alcentr"><i>Критерій<br>непридатний</i></td><td class="alcentr"><i>Критерій<br>непридатний</i></td><td class="alcentr"><i>Критерій<br>непридатний</i></td></tr>
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
  <li>Збірка генома утримує принаймні один запис бази даних RefSeq, у назві якої є фраза «map unlocalized», і цей запис утримує послідовность гена 16S рРНК або його частину.</li>
</ul>
<br>
<p>Програми, за допомогую яких RiboGrove була створена, знаходиться у наступному ҐітХаб-репозиторії: <a href="https://github.com/masikol/ribogrove-tools">ribogrove-tools</a>.</p>
<hr>
{% endif %}
<div id="downloads" class="pad-anchor"></div>
<h2>Завантаження</h2>
<div id="current-release" class="pad-anchor"></div>
{% if not archive %}
<h3>Актуальний реліз RiboGrove — {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
{% else %}
<h3>Реліз RiboGrove {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
{% endif %}
<p>Реліз заснований на базі даних RefSeq {{ refseq_release }}.</p>
<ul>
<li>Fasta-файл повнорозмірних послідовностей генів 16S рРНК. <a href="ribogrove_releases/{{ ribogrove_release_number }}/ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Завантажити („ґзіпʼнутий“ fasta-файл, {{ final_fasta_fsize_fmt }} Мб)</a></li>

<li>Метадані. <a href="ribogrove_releases/{{ ribogrove_release_number }}/metadata_ribogrove_{{ ribogrove_release_number }}.zip">Завантажити (zip-архів {{ metadata_fsize_fmt }} Мб)</a><br>
<details>
  <summary><small>Яка саме інформація утримується у метаданих?</small></summary>
  <small>
    <p>Метадані складаються з наступних файлів:</p>
    <ol>
    <li>
      <span class="samp-highl">discarded_sequences.fasta.gz</span><br>
      Це fasta-файл, який утримує послідовності, які були у геномах-джерелах, були при цьому анотувані як гени 16S рРНК, але які були „забракувані“ (discarded) через їхню неповноту, присутність внутрішньогенних повторів і т.д. Отже, ці послідовності не включені до RiboGrove.
    </li>
    <li>
      <span class="samp-highl">source_RefSeq_genomes.tsv</span><br>
      Це TSV-файл, який утримує інформацію про те, які геноми выкористовувалися для створення RiboGrove.
    </li>
    <li>
      <span class="samp-highl">gene_seqs_statistics.tsv</span>, <span class="samp-highl">discarded_gene_seqs_statistics.tsv</span><br>
      Це TSV-файли, які утримують інформацію пра нуклеотидний склад, розмір, геномну і таксономічну приналежність послідовностей генів. Перший файл із оцих двох описує фінальні послідовності, включені до RiboGrove, а другий — „забракувані“ послідовності.
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
      <span class="samp-highl">entropy_summary.tsv</span><br>
      Це TSV-файл, який утримує інформацію про внутрішньогеномну мінливість генів 16S рРНК. Внутрішньогеномну мінливість розраховували тільки для геномів категорії&nbsp;1, у яких є принаймні два гени 16S рРНК. Для оцінки мінливості выкористовували ентропію Шенона: вирівнювали послідовності генів кожного геному за допомогою програми <a href="https://drive5.com/muscle5/">MUSCLE</a>, і після розраховували ентропію для кожної колонки множинного вирівнювання.
    </li>
    <li>
      <span class="samp-highl">primer_pair_genomic_coverage.tsv</span><br>
      Це TSV-файл, який утримує геномне покриття пар ПЛР-праймерів до різних V-регіонів генів 16S рРНК. Наприклад, для родини <i>Enterobacteriaceae</i> геномне покриття пари праймерів — це відсоток геномів <i>Enterobacteriaceae</i>, які утримують принаймні один ген 16S рРНК, на матриці якого (теоретично) може сформуватися ПЛР-продукт за допомогою цієї пари праймерів.
    </li>
  </ol></small>
</details>
</li>
</ul>
<p>Файл fasta стиснутий за допомогою програми gzip, а файли метаданих стиснуті у zip-архів. Щоби розпакувати їх, користувачі Linux та Mac OS можуть скористатися саме програмамі gzip та zip; оці програми звичайно установлені разом із операційною системою. Користувачам Windows може придатися безкоштовна програма-архіватор <a href="https://www.7-zip.org/">7-Zip</a>.</p>
<div id="release-archive" class="pad-anchor"></div>
<h3>Архів релізів RiboGrove</h3>
<p>Всі попередні релізи можна знайти в <a href="archive/ribogrove_archive_en.html"> архіві релізів RiboGrove</a>. На жаль, він ведеться тільки англійською.</p>
<div id="release-notes" class="pad-anchor"></div>
<h3>Коментарі до релізу</h3>
<p><i>Цей реліз нічим принципово не відрізнюється від попереднього.</i></p>
<p>Коментарі до всіх релізів RiboGrove можна знайти на <a href="ribogrove_release_notes_en.html">сторінці коментарів до релізів</a>.</p>
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
<tr><th class="numcol">Кількість копій <sup>*</sup></th><th class="alncenter" colspan="2">Бактерії</th><th class="alncenter" colspan="2">Археї</th></tr>
<tr><td></td><td>Кількість видів</td><td>Відсоток видів, %</td><td>Кількість видів</td><td>Відсоток видів, %</td></tr>
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
<p><sup>*</sup> Медіанна кількість копій у межах виду.</p>

<div id="top-longest" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 найдовших генів 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Організм</th><th class="numcol">Довжина, п.о.</th><th class="alnleft">SeqID <sup>*</sup></th><th class="alnleft">Код доступу<br>геномної збірки</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерії</td>
</tr>
{% for _, row in ribogrove_top_longest_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Археї</td>
</tr>
{% for _, row in ribogrove_top_longest_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>
<p><sup>*</sup> SeqID — це ідентифікатор послідовності RiboGrove.</p>

<div id="top-shortest" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 найкоротших генів 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Організм</th><th class="numcol">Довжина, п.о.</th><th class="alnleft">SeqID <sup>*</sup></th><th class="alnleft">Код доступу<br>геномної збірки</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерії</td>
</tr>
{% for _, row in ribogrove_top_shortest_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Археї</td>
</tr>
{% for _, row in ribogrove_top_shortest_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>
<p><sup>*</sup> SeqID — це ідентифікатор послідовності RiboGrove.</p>

<div id="top-copy-num" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 геномів за найбільшею кілкістю копій генів 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Організм</th><th class="numcol">Кількість копій</th><th class="alnleft">Код доступу<br>геномної збірки</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерії</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Археї</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<div id="top-var" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 геномів за найбільшею мінливістю генів 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Організм</th><th class="numcol">Сума ентропії<sup> *</sup>, біти</th><th class="numcol">Середня ентропія<sup> *</sup>, біти</th><th class="numcol">Кількість мінливих позіцій</th><th class="numcol">Кількість копій генів</th><th class="alnleft">Код доступу<br>геномної збірки</th></tr>
<tr>
<td colspan="6" class="subhead">Бактерії</td>
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
<td colspan="6" class="subhead">Археї</td>
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

<p><sup>*</sup> Ентропія тут — це <a href="https://uk.wikipedia.org/wiki/Інформаційна_ентропія">ентропія Шеннона</a>, яку розрахувовали для кожної колонки множинного вирівнювання повнорозмірних послідовностей генів 16S рРНК кожного геному. І далі ентропію підсумовували (колонка «Сума ентропії») і усереднювали (колонка «Середня ентропія»).</p>
<p><sup>**</sup> <i>Halomicrobium</i> sp. ZPS1 — показовий випадок. У цьому геномі два гени 16S рРНК, тому ентропія є рівною кількості нуклеотидів, які не співпадають у їхніх послідовностях. Коли підрахувати відсоток ідентичності (скорочено ВІ) послідовностей оцих двох генів, отримуємо 90,70%! Нагадуємо, що за поріг ВІ послідовностей генів 16S рРНК для віднесення організмів до відрізных <em>родів</em> звычайно береться (хоча і довільно) значення 95%.</p>

<div id="primers-coverages" class="pad-anchor"></div>
<table class="sum-table"><caption>Спектр диї<sup>*</sup> пар праймерів до різних V-регіонів бактеріальних генів 16S рРНК</caption>
<tbody class="primer-cov-tbody">
<tr>
  <th class="alnleft" rowspan="2">Тип<br />бактерий</th>
  <th class="numcol" rowspan="2">Кількість<br />геномів</th>
  <th class="numcol">Повний<br />ген</th>
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

<p><sup>*</sup> Спектр диї пары праймерів — це видсоток геномів, які мають принаймні один ген 16S рРНК, який можливо заампліфікувати за допомогою ПЦР і з застосуванням цієї пари праймерів. Деталі описані у нашому <a href="https://masikol.github.io/latest_ribogrove_release_ua.html#citing-ribogrove">артикулі про RiboGrove</a>.</p>
<p>Є більш детальна таблиця — у файлі <span class="samp">primer_pair_genomic_coverage.tsv</span> у метаданих <a href="#downloads">metadata</a>. Та таблиця утримує покриття не тільки для типів бактерий, але і для кожного класу, ряду, родини, роду та виду. Крім того, та таблиця утримує покриття для пари праймерів 1115F–1492R (регіон V7–V9). У цієї таблиці нема цієї пари, щоб таблиця не була занадто великою.</p>

<table class="sum-table"><caption>Праймери, спектр диї яких оцінювався</caption>
<tbody>
  <tr><th>Назва праймера</th><th>Послідовність</th><th>Посилання</th></tr>
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

<p><sup>*</sup> Праймери 341F та 785R застосовуються у <a href="https://support.illumina.com/downloads/16s_metagenomic_sequencing_library_preparation.html">протоколі</a> підготовки бібліотеки для секвенування ампліконів регіонів V3–V4 генів 16S рРНК на апараті Illumina MiSeq.</p>

<hr>
<div id="searching-data" class="pad-anchor"></div>
<h2>Пошук даних у RiboGrove</h2>
<p>RiboGrove — це дуже мінімалістична база даних: вона складається з двох файлів fasta із додатковими метаданими. Отже, розширених засобів пошуку немає. Ми визнаємо цю проблему і залишаємо нижче делілька рекомендацій щодо пошуку даних у RiboGrove.</p>
<div id="header-format" class="pad-anchor"></div>
<h3>Формат заголовку</h3>
<p>Записи послідовностей RiboGrove у форматі fasta мають наступний формат заголовку:</p>
<p class="samp-highl samp-vwide" style="font-size:14px">&gt;GCF_000978375.1:NZ_CP009686.1:8908-10459:plus ;d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__cereus; category:1</p>
<p>Основні частини заголовку відокремлені один від одного пробілами. В заголовку повинні бути три такі частини:</p>
<ol>
<li>Ідентифікатор послідовності (SeqID): <span class="samp-highl">GCF_000978375.1:NZ_CP009686.1:8908-10459:plus</span>. SeqID, у свою чергу, складається з чотирьох частинок, відокремлених двокрапками:
<ol>
<li>Код доступу збірки геному, від якого походить послідовність гену: <span class="samp-highl">GCF_000978375.1</span>.</li>
<li>Код доступу послідовності RefSeq, від якої походить послідовність гену: <span class="samp-highl">NZ_CP009686.1</span>.</li>
<li>Координати гену у цієї геномної послідовності RefSeq: <span class="samp-highl">8908-10459</span> (нумерація йде з 1-ці, інтервали включають ліві та праві межі).</li>
<li>Ланцюг геномної послідовності, на якій ген знаходиться: <span class="samp-highl">plus</span> (або <span class="samp-highl">minus</span>).</li>
</ol>
</li>
<li>Рядок таксономії, а саме назви домену (<span class="samp-highl">Bacteria</span>), типу (<span class="samp-highl">Firmicutes</span>), класу (<span class="samp-highl">Bacilli</span>), ряду (<span class="samp-highl">Bacillales</span>), родини (<span class="samp-highl">Bacillaceae</span>) та роду (<span class="samp-highl">Bacillus</span>), а також видовий епітет (<span class="samp-highl">cereus</span>).<br>
Перед кожній назвою пишеться префікс, який позначує ранг: <span class="samp-highl">d__</span> для домену, <span class="samp-highl">p__</span> для типу, <span class="samp-highl">c__</span> для класу, <span class="samp-highl">o__</span> для ряду, <span class="samp-highl">f__</span> для родини, <span class="samp-highl">g__</span> для роду і <span class="samp-highl">s__</span> для видового епітету. У префіксах після літеры пишуться <strong>два</strong> знаки нижнього підкресливання.<br>
Таксономічні назви відокремлені та оточені крапками з комою (<span class="samp-highl">;</span>).</li>
<li>Категорія геному, із якого послідовність гену походить: (<span class="samp-highl">category:1</span>).</li>
</ol>
<h3>Вибір послідовностей</h3>
<p>Відбирати ті чи інші послідовності із файлів fasta можна за допомогою програми <strong>Seqkit</strong> (<a href="https://github.com/shenwei356/seqkit">ҐітХаб-репозиторій</a>, <a href="https://bioinf.shenwei.me/seqkit/">документація</a>). Це безкоштовна, кросплатформенна, багатофункціональна та доволі швидка програма, яка може обробляти як „ґзіпʼнуті“, так і прості нестиснуті файли fasta. Для відбіру послідовностей зручно выкористовувати програми <span class="samp">seqkit grep</span> та <span class="samp">seqkit seq</span>.</p>
<div id="subset-by-header" class="pad-anchor"></div>
<h4>Вибір послідовностей по заголовку</h4>
<p>Нехай ми робимо вибірку з файла <span class="samp">ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</span>. Корисними можуть статися наступні приклади команд користання програмою <span class="samp">seqkit grep</span>:</p>
<p><strong>Приклад 1</strong>. Вибрати одну послідовність по SeqID.</p>
<p class="samp-highl samp-vwide">seqkit grep -p "GCF_000978375.1:NZ_CP009686.1:8908-10459:plus" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опція <span class="samp">-p</span> задає рядок, який програма буде шукати у fasta-заголовках (на самій справі, тільки в тій частині заголовків, де запісаний SeqID).</p>
<p><strong>Приклад 2</strong>. Вибрати всі послідовності генів однієї геномної послідовності RefSeq по коду доступу <span class="samp">NZ_CP009686.1</span>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ":NZ_CP009686.1:" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Тут потрібні ще дві опціі: <span class="samp">-n</span> та <span class="samp">-r</span>. Перша вказує програмі шукати підрядки в заголовках цілком, а не тільки в ідентифікаторах (SeqID). Друга опція вказує програмі, що їй треба шукати не тільки ті заголовки, які цілком співпадають із шуканим рядком, а також і ті, які утримують шуканий рядок як свою частину.</p>
<p class="samp-comment">Щоб забезпечити специфичність пошуку, оточюйте код доступу двокрапками (<span class="samp-highl">:</span>).</p>
<p><strong>Приклад 3</strong>. Вибрати всі послідовності генів одного геному (код доступу геномної збірки <a href="https://ncbi.nlm.nih.gov/datasets/genome/GCF_019357495.1">GCF_019357495.1</a>).</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "GCF_019357495.1:" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Щоб забезпечити специфичність пошуку, пишіть двокрапку (<span class="samp-highl">:</span>) після коду доступу збірки.</p>
<p><strong>Приклад 4</strong>. Вибрати всі послідовності актинобактерий.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ";p__Actinobacteria;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Щоб забезпечити специфичність пошуку, оточюйте таксономію крапками з комою (<span class="samp-highl">;</span>).</p>
<p><strong>Приклад 5</strong>. Вибрати всі послідовності з геномів категорії 1.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "category:1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Приклад 6</strong>. Вибрати всі послідовності крім тих, які належать до типу <i>Firmicutes</i>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nvrp ";p__Firmicutes;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Зверніть увагу на опцію <span class="samp">-v</span> у рядку з опціями <span class="samp">-nvrp</span>. Оця опція інвертує поведінку програми, і тепер вона знайде всі послідовності, заголовки яких <em>не</em> утримують рядку «<span class="samp">;p__Firmicutes;</span>».</p>
<div id="subset-by-len" class="pad-anchor"></div>
<h4>Вибір послідовностей по довжині</h4>
<p>Відбірати послідовності по довжині можна за допомогою програми <span class="samp">seqkit seq</span>.</p>
<p><strong>Приклад 1</strong>. Вибрати всі послідовності довше ніж 1600 п.о.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1601 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опція <span class="samp">-m</span> задає мінімальну довжину послідовностей, які програма подасть на вихід.</p>
<p><strong>Приклад 2</strong>. Вибрати всі послідовності коротше ніж 1500 п.о.</p>
<p class="samp-highl samp-vwide">seqkit seq -M 1499 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опція <span class="samp">-M</span> задає максимальну довжину послідовностей, які програма подасть на вихід.</p>
<p><strong>Приклад 3</strong>. Вибрати всі послідовності, довжина яких знаходиться в інтервалі [1500, 1600] п.о.</p>
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
<p><strong>Приклад 3</strong>. Вибрати всі коди доступу геномних послідовностей з БД RefSeq.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f2 -d':' | sort | uniq</p>
<p class="samp-comment">Така команда спрацює, тільки якщо на компʼютері установлені утиліти <span class="samp">cut</span>, <span class="samp">sort</span>, та <span class="samp">uniq</span> (на Linux та Mac OS оці програми звичайно встановлені разом із сістемою).</p>

<p><strong>Приклад 4</strong>. Вибрати всі ідентифікатори геномних збірок.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f1 -d':' | sort | uniq</p>
<p class="samp-comment">Така команда спрацює, тільки якщо на компʼютері установлені утиліти <span class="samp">cut</span>, <span class="samp">sed</span>, <span class="samp">sort</span> і <span class="samp">uniq</span> (на Linux та Mac OS оці програми звичайно встановлені разом із сістемою).</p>
<p><strong>Приклад 5</strong>. Вибрати всі назви типів організмів.</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | grep -Eo ';p__[^;]+' | sed -E 's/;|p__//g' | sort | uniq</p>
<p class="samp-comment">Така команда спрацює, тільки якщо на компʼютері установлені утиліти <span class="samp">grep</span>, <span class="samp">sed</span>, <span class="samp">sort</span> і <span class="samp">uniq</span> (на Linux та Mac OS оці програми звичайно встановлені разом із сістемою).</p>
<hr>
<div id="contacts" class="pad-anchor"></div>
<h2>Контакти</h2>
<p>Якщо ви маєте питання щодо RiboGrove, будь ласка, звертайтеся до Максима Сиколенка на ел. адресу sikolenko<img class="sabaka" src="img/pes2.gif" alt="[ at ]" align="bottom">mbio.bas-net.by або maximdeynonih<img class="sabaka" src="img/pes2.gif" alt="[ at ]" align="bottom">gmail.com.</p>
{% if not archive %}
<hr>
<div id="citing-ribogrove" class="pad-anchor"></div>
<h2>Цитування RiboGrove</h2>
<p>Якщо база даних RiboGrove сталася корисною для ваших досліджень, посилайтеся, будь ласка, на наступну публікацію:</p>
<p>Maxim A. Sikolenko, Leonid N. Valentovich. “RiboGrove: a database of full-length prokaryotic 16S rRNA genes derived from completely assembled genomes” // Research in Microbiology, Volume 173, Issue 4, May 2022, 103936.<br> (ЦІО: <a href="https://doi.org/10.1016/j.resmic.2022.103936">10.1016/j.resmic.2022.103936</a>).</p>
<hr>
<div id="faq" class="pad-anchor"></div>
<h2>Питання, які люди задають про RiboGrove</h2>
<div><details>
  <summary><b>1. Як зробити файл з таксономією, з яким могли б працювати програми з пакета QIIME2?</b></summary>
  <p>Для цього можете скористатися скриптом <a href="https://github.com/masikol/ribogrove-tools/blob/main/data_transformation_scripts/make_qiime_taxonomy_file.py"><span class="samp">make_qiime_taxonomy_file.py</span></a>: він трансформує файл <span class="samp">metadata/taxonomy.tsv</span> в файл, з яким змогуть працювати програми QIIME2. Довідку-керівництво для цього скрипта можете знайти у відповідним <a href="https://github.com/masikol/ribogrove-tools/blob/main/data_transformation_scripts/README.md">файле README</a>.</p>
</details></div><br>
<div><details>
  <summary><b>2. Як зберегти у файлі дані, відібрані за допомогою Seqkit?</b></summary>
  <p>У відповідної дискусії на форумі вже дали декілька корисних відповідей та порад: <a href="https://bioinformatics.stackexchange.com/questions/20915/how-do-i-save-selected-sequences-in-seqkit-to-a-file">https://bioinformatics.stackexchange.com/questions/20915/how-do-i-save-selected-sequences-in-seqkit-to-a-file</a>.</p>
</details></div><br>
<div><details>
  <summary><b>3. Як знайти потрібні записи у fasta-файлі по послідовності за допомогою Seqkit?</b></summary>
  <p>У відповідної дискусії на форумі вже дали декілька корисних відповідей та порад: <a href="https://www.biostars.org/p/9561418">https://www.biostars.org/p/9561418</a>.</p>
</details></div>
{% endif %}
<br>
<hr>
<p>RiboGrove, {{ ribogrove_release_date }}</p>
</div>
</body>
</html>
