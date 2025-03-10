<h1 style="display:none">RiboGrove</h1>
<img src="/cager/img/RiboGrove_logo.png"><br><br>
<h2>Содержание</h2>
<ul>
{% if not archive %}
<li><a href="#overview">Что такое RiboGrove</a>
<ul class="ribogrove-nested-list">
<li><a href="#other-databases">RiboGrove и другие базы данных 16S рРНК</a></li>
<li><a href="#genome-categories">Категории геномов</a></li>
</ul>
</li>
{% endif %}
<li><a href="#downloads">Файлы</a>
<ul class="ribogrove-nested-list">
  {% if not archive %}
  <li><a href="#current-release">Актуальный выпуск RiboGrove — {{ ribogrove_release_number }}</li>
  {% else %}
  <li><a href="#current-release">Выпуск RiboGrove {{ ribogrove_release_number }}</li>
  {% endif %}
  <li><a href="#release-archive">Архив выпусков RiboGrove</li>
  <li><a href="#release-notes">Комментарии к выпуску</li>
</ul>
</li>
<li><a href="#summary">Статистическое описание</a>
<ul class="ribogrove-nested-list">
<li><a href="#db-size">Размер базы данных RiboGrove</a></li>
<li><a href="#gene-lens">Длина генов 16S рРНК</a></li>
<li><a href="#copy-number">Число копий генов 16S рРНК</a></li>
<li><a href="#top-longest">Топ-10 самых длинных генов 16S рРНК</a></li>
<li><a href="#top-shortest">Топ-10 самых коротких генов 16S рРНК</a></li>
<li><a href="#top-copy-num">Tоп-10 геномов с наибольшим числом копий генов 16S рРНК</a></li>
<li><a href="#top-var">Tоп-10 геномов с наибольшей внутригеномной изменчивостью генов 16S рРНК</a></li>
<li><a href="#primers-coverages">Спектр действия пар праймеров к разным V-участкам бактериальных генов 16S рРНК</a></li>
</ul>
</li>
<li><a href="#searching-data">Поиск данных в базе данных RiboGrove</a>
<ul class="ribogrove-nested-list">
<li><a href="#header-format">Формат заголовка</a></li>
<li><a href="#subset-by-header">Выборка последовательностей по заголовку</a></li>
<li><a href="#subset-by-len">Выборка последовательностей по длине</a></li>
<li><a href="#select-head">Извлечение информации из заголовков последовательностей</a></li>
</ul>
</li>
{% if not archive %}
<li><a href="#contacts">Контакты</a></li>
<li><a href="#citing-ribogrove">Цитирование RiboGrove</a></li>
<li><a href="#faq">Вопросы, которые возникают у людей про RiboGrove</a></li>
{% endif %}
</ul>
<hr>
{% if not archive %}
<div id="overview" class="pad-anchor"></div>
<h2>Что такое RiboGrove</h2>
<p>RiboGrove — это база данных последовательностей генов 16S рРНК бактерий и архей.</p>
<p>RiboGrove основывается на базе данных <a href="https://ncbi.nlm.nih.gov/refseq/">RefSeq</a>. RiboGrove содержит только полноразмерные последовательности генов 16S рРНК, а сами последовательности взяты из полностью собранных геномов прокариот, депонированных в RefSeq. Поэтому мы постулируем высокую надёжность последовательностей, содержащихся в RiboGrove.</p>
<div id="other-databases" class="pad-anchor"></div>
<h3>RiboGrove и другие базы данных 16S рРНК</h3>
<p>Таблица ниже показывает, чем (качественно) RiboGrove отличается от схожих баз данных последовательностей РНК, а именно <a href="https://rrndb.umms.med.umich.edu/"><i>rrn</i>DB</a>, <a href="https://www.arb-silva.de/">Silva</a>, <a href="http://rdp.cme.msu.edu/">RDP</a> и <a href="https://greengenes.secondgenome.com/">Greengenes</a>. Говоря коротко, RiboGrove уступает аналогам по количеству последовательностей и их разнообразию, но превосходит по достоверности последовательностей.</p>
<table class="sum-table">
<tbody>
  <tr class="sumtab-row"><th></th><th class="alcentr">RiboGrove</th><th class="alcentr"><i>rrn</i>DB</th><th class="alcentr">Silva</th><th class="alcentr">RDP</th><th class="alcentr">Greengenes</th></tr>
  <tr class="sumtab-row"><td class="feature-col">Представленные организмы</td><td class="alcentr">Бактерии<br>Археи</td><td class="alcentr">Бактерии<br>Археи</td><td class="alcentr">Бактерии<br>Археи<br>Эукариоты</td><td class="alcentr">Бактерии<br>Археи<br>Эукариоты</td><td class="alcentr">Бактерии<br>Археи</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Представленные субъединицы рибосом</td><td class="alcentr">Малая</td><td class="alcentr">Малая</td><td class="alcentr">Большая<br>Малая</td><td class="alcentr">Большая<br>Малая</td><td class="alcentr">Малая</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Содержит последовательности собранных геномов</td><td class="alcentr">Да</td><td class="alcentr">Да</td><td class="alcentr">Да</td><td class="alcentr">Да</td><td class="alcentr">Да</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Содержит последовательности ампликонов</td><td class="alcentr">Нет</td><td class="alcentr">Нет</td><td class="alcentr">Да</td><td class="alcentr">Да</td><td class="alcentr">Да</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Содержит неполные последовательности генов</td><td class="alcentr">Нет</td><td class="alcentr">Да</td><td class="alcentr">Да</td><td class="alcentr">Да</td><td class="alcentr">Да</td></tr>
  <tr class="sumtab-row"><td class="feature-col">Различает категории геномов</td><td class="alcentr">Да</td><td class="alcentr">Нет</td><td class="alcentr"><i>Критерий<br>неприменим</i></td><td class="alcentr"><i>Критерий<br>неприменим</i></td><td class="alcentr"><i>Критерий<br>неприменим</i></td></tr>
</tbody>
</table>
<div id="genome-categories" class="pad-anchor"></div>
<h3>Категории геномов</h3>
<p>Все геномы, из которых брались данные для создания RiboGrove, были разделены на три категории по их достоверности:</p>
<ol style="list-style:none;">
  <li><strong>Категория 1</strong> (наибольшая достоверность). Геномы, которые не демонстрируют признаков низкокачественной сборки и которые были секвенированы с помощью технологии PacBio либо комбинации (Oxford&nbsp;Nanopore&nbsp;+&nbsp;Illumina).</li>
  <li><strong>Категория 2.</strong> Геномы, которые не демонстрируют признаков низкокачественной сборки, секвенированные с помощью любой другой технологии (либо те, для которых технология секвенирования не была указана).</li>
  <li><strong>Категория 3</strong> (наименьшая достоверность). Геномы, которые демонстрируют хотя бы один признак низкокачественной сборки.</li>
</ol>
<p>Признаки низкокачественной сборки:</p>
<ul>
  <li>Геном содержит хотя бы одно вырожденное основание в последовательностях его генов 16S рРНК.</li>
  <li>Сборка генома содержит хотя бы одну запись базы данных RefSeq, в названии которой есть фраза «map unlocalized», и эта запись содержит последовательность гена 16S рРНК либо его часть.</li>
</ul>
<br>
<p>Программы, с помощью которых RiboGrove была создана, доступны в репозитории на ГитХабе: <a href="https://github.com/masikol/ribogrove-tools">ribogrove-tools</a>.</p>
<hr>
{% endif %}
<div id="downloads" class="pad-anchor"></div>
<h2>Файлы</h2>
<div id="current-release" class="pad-anchor"></div>
{% if not archive %}
<h3>Актуальный выпуск RiboGrove — {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
{% else %}
<h3>Выпуск RiboGrove {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
{% endif %}
<p>Выпуск RiboGrove {{ ribogrove_release_number }} основан на данных, полученных из базы данных RefSeq версии {{ refseq_release }}.</p>
<ul>
<li>Fasta-файл полноразмерных последовательностей генов 16S рРНК. <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Скачать (gzipʼнутый fasta-файл, {{ final_fasta_fsize_fmt }} Мб)</a></li>
<li>Метаданные. <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/metadata_ribogrove_{{ ribogrove_release_number }}.zip">Скачать (zip-архив {{ metadata_fsize_fmt }} Мб)</a><br>
<details>
  <summary><small>Какая именно информация содержится в метаданных?</small></summary>
  <small>
    <p>Метаданные содержат следующие файлы:</p>
    <ol>
    <li>
      <span class="samp-highl">discarded_sequences.fasta.gz</span><br>
      Это fasta-файл, который содержит последовательности, которые присутствовали в исходных геномах, были пры этом аннотированы как гены 16S рРНК, но которые были „забракованы“ (discarded) из-за их неполноты, наличия внутригенных повторов и т.д. Таким образом, эти последовательности не включены в RiboGrove.
    </li>
    <li>
      <span class="samp-highl">source_RefSeq_genomes.tsv</span><br>
      Это TSV-файл, который содержит информацию про то, какие геномы использовались для создания RiboGrove.
    </li>
    <li>
      <span class="samp-highl">gene_seqs_base_counts.tsv</span>, <span class="samp-highl">discarded_gene_seqs_base_counts.tsv</span><br>
      Это TSV-файлы, которые содержат информацию про нуклеотидный состав и размер последовательностей генов. Первый файл их этих двух описывает финальные последовательности, включенные в RiboGrove, а второй — „забракованые“ последовательности.
    </li>
    <li>
      <span class="samp-highl">categories.tsv</span><br>
      Это TSV-файл, который содержит информацию про категорию каждого генома, который был использован для создания RiboGrove. Также в этом файле есть информация про технологии, которые использовались для секвенирования этих геномов.
    </li>
    <li>
      <span class="samp-highl">taxonomy.tsv</span><br>
      Это TSV-файл, который содержит таксономическую принадлежность каждого генома и гена.
    </li>
    <li>
      <span class="samp-highl">intragenic_repeats.tsv</span><br>
      Это TSV-файл, который содержит информацию про внутригенные повторы в последовательностях генов, обнаруженные с помощью программы <a href="https://github.com/deprekate/RepeatFinder">RepeatFinder</a>.
    </li>
    <li>
      <span class="samp-highl">entropy_summary.tsv</span><br>
      Это TSV-файл, который содержит информацию про внутригеномную изменчивость генов 16S рРНК. Внутригеномную изменчивость вычисляли только для геномов категории&nbsp;1, в которых есть хотя бы два гена 16S рРНК. Для оценки изменчивости использовали энтропию Шеннона: выравнивали последовательности генов каждого генома с помощью программы <a href="https://mafft.cbrc.jp/alignment/software/">MAFFT</a>, и после этого вычисляли энтропию для каждого столбца множественного выравнивания.
    </li>
    <li>
      <span class="samp-highl">primer_pair_genomic_coverage.tsv</span><br>
      Это TSV-файл, который содержит геномное покрытие пар ПЦР-праймеров к разным V-участкам генов 16S рРНК. Например, для семейства <i>Enterobacteriaceae</i> геномное покрытие пары праймеров — это процент геномов <i>Enterobacteriaceae</i>, которые содержат как минимум один ген 16S рРНК, на матрице которого (теоритично) может сформироваться ПЦР-продукт с помощью этой пары праймеров.
    </li>
  </ol></small>
</details>
</li>
</ul>
<p>Fasta-файл сжат с помощью программы gzip, а метаданные — с помощью программы zip. Чтобы разархивировать их, пользователи Linux и Mac OS могут воспользоваться программами gzip и zip, соответственно; эти программы должны быть предустановлены. Пользователи же Windows могут воспользоваться бесплатной программой <a href="https://www.7-zip.org/">7-Zip</a>.</p>
<div id="release-archive" class="pad-anchor"></div>
<h3>Архив выпусков RiboGrove</h3>
<p>Вы можете найти все выпуски RiboGrove в <a href="/cager/ru/content/60-ribogrove-archive">архиве выпусков RiboGrove</a>.</p>
<div id="release-notes" class="pad-anchor"></div>
<h3>Комментарии к выпуску</h3>
<p><i>Нет принципиальных отличий от предыдущего выпуска.</i></p>
<p>Комментарии ко всем выпускам RiboGrove можно найти на <a href="/cager/ru/content/66-ribogrove-release-notes">странице комментариев к выпускам</a>.</p>
<hr>
<div id="summary" class="pad-anchor"></div>
<h2>Статистическое описание</h2>
<div id="db-size" class="pad-anchor"></div>
<table class="sum-table">
<caption>Размер базы данных RiboGrove</caption>
<tbody>
<tr><th></th><th class="numcol">Бактерии</th><th class="numcol">Археи</th><th class="numcol">Всего</th></tr>
<tr class="sumtab-row">
<td>Количество последовательностей генов</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['gene_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Количество уникальных последовательностей генов</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['uniq_gene_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Количество видов</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['species_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Количество геномов</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Количество геномов 1-й категории</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat1_genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Количество геномов 2-й категории</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat2_genome_num']['Total'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Количество геномов 3-й категории</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Archaea'] }}</td>
<td class="numcol">{{ ribogrove_size_dict['cat3_genome_num']['Total'] }}</td>
</tr>
</tbody>
</table>
<div id="gene-lens" class="pad-anchor"></div>
<table class="sum-table">
<caption>Длина генов 16S рРНК</caption>
<tbody>
<tr><th></th><th class="numcol">Бактерии</th><th class="numcol">Археи</th></tr>
<tr class="sumtab-row">
<td>Минимум, п.н.</td>
<td class="numcol">{{ ribogrove_len_dict['min']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['min']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>25-й перцентиль, п.н. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['25perc']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['25perc']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Медиана, п.н. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['median']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['median']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>75-й перцентиль, п.н. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['75perc']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['75perc']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Средняя длина, п.н. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['mean']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['mean']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Мода, п.н. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['modes']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['modes']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Максимальная длина, п.н.</td>
<td class="numcol">{{ ribogrove_len_dict['max']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['max']['Archaea'] }}</td>
</tr>
<tr class="sumtab-row">
<td>Среднеквадратическое отклонение, п.н. <sup>*</sup></td>
<td class="numcol">{{ ribogrove_len_dict['std']['Bacteria'] }}</td>
<td class="numcol">{{ ribogrove_len_dict['std']['Archaea'] }}</td>
</tr>
</tbody>
</table>
<p><sup>*</sup> Метрики, помеченные звёздочкой, были рассчитаны с предварительной нормализацией данных, а именно: медианные значения длины внутри каждого вида использовались для расчёта этих метрик.</p>
<div id="copy-number" class="pad-anchor"></div>
<table class="sum-table">
<caption>Число копий генов 16S рРНК</caption>
<tbody>
<tr><th class="numcol">Число копий <sup>*</sup></th><th class="alncenter" colspan="2">Бактерии</th><th class="alncenter" colspan="2">Археи</th></tr>
<tr><td></td><td>Количество видов</td><td>Процент видов, %</td><td>Количество видов</td><td>Процент видов, %</td></tr>
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
<p><sup>*</sup> Медианное число копий в пределах вида.</p>

<div id="top-longest" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 самых длинных генов 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Организм</th><th class="numcol">Длина гена, п.н.</th><th class="alnleft">SeqID <sup>*</sup></th><th class="alnleft">Код доступа<br>геномной сборки</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерии</td>
</tr>
{% for _, row in ribogrove_top_longest_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name'], row['asm_acc']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Археи</td>
</tr>
{% for _, row in ribogrove_top_longest_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name'], row['asm_acc']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>
<p><sup>*</sup> SeqID — это идентификатор последовательности RiboGrove.</p>

<div id="top-shortest" class="pad-anchor"></div>
<table class="sum-table">
<caption>Топ-10 самых коротких генов 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Организм</th><th class="numcol">Длина гена, п.н.</th><th class="alnleft">SeqID <sup>*</sup></th><th class="alnleft">Код доступа<br>геномной сборки</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерии</td>
</tr>
{% for _, row in ribogrove_top_shortest_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name'], row['asm_acc']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Археи</td>
</tr>
{% for _, row in ribogrove_top_shortest_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name'], row['asm_acc']) }}</td>
<td class="numcol">{{ row['len'] }}</td>
<td class="seqid-td">{{ '<br>'.join(row['seqID']) }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>
<p><sup>*</sup> SeqID — это идентификатор последовательности RiboGrove.</p>

<div id="top-copy-num" class="pad-anchor"></div>
<table class="sum-table">
<caption>Tоп-10 геномов с наибольшим числом копий генов 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Организм</th><th class="numcol">Число копий</th><th class="alnleft">Код доступа<br>геномной сборки</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерии</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name'], row['asm_acc']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Археи</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name'], row['asm_acc']) }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://ncbi.nlm.nih.gov/datasets/genome/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<div id="top-var" class="pad-anchor"></div>
<table class="sum-table">
<caption>Tоп-10 геномов с наибольшей внутригеномной изменчивостью генов 16S рРНК</caption>
<tbody>
<tr><th class="alnleft">Организм</th><th class="numcol">Сумма энтропии<sup> *</sup>, биты</th><th class="numcol">Средняя энтропия<sup> *</sup>, биты</th><th class="numcol">Количество изменчивых позиций</th><th class="numcol">Число копий генов</th><th class="alnleft">Код доступа<br>геномной сборки</th></tr>

<tr>
<td colspan="6" class="subhead">Бактерии</td>
</tr>
{% for _, row in ribogrove_top_intragenomic_var_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name'], row['asm_acc']) }}</td>
<td class="numcol">{{ row['sum_entropy'] }}</td>
<td class="numcol">{{ row['mean_entropy'] }}</td>
<td class="numcol">{{ row['num_var_cols'] }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://www.ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
<tr>
<td colspan="6" class="subhead">Археи</td>
</tr>
{% for _, row in ribogrove_top_intragenomic_var_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td>{{ retrieve_strain_name(row['strain_name'], row['asm_acc']) }}</td>
<td class="numcol">{{ row['sum_entropy'] }}</td>
<td class="numcol">{{ row['mean_entropy'] }}</td>
<td class="numcol">{{ row['num_var_cols'] }}</td>
<td class="numcol">{{ row['copy_number'] }}</td>
<td><a href="https://www.ncbi.nlm.nih.gov/assembly/{{ row['asm_acc'] }}">{{ row['asm_acc'] }}</a></td>
</tr>
{% endfor %}
</tbody>
</table>

<p><sup>*</sup> Энтропия в данном случае — это <a href="https://ru.wikipedia.org/wiki/Информационная_энтропия">энтропия Шеннона</a>, вычисленная для каждого столбца множественного выравнивания полноразмерных последовательностей генов 16S рРНК индивидуального генома. Вычисленная энтропия далее суммировалась (столбец «Сумма энтропии») или усреднялась (столбец «Средняя энтропия»).</p>
<p><sup>**</sup> <i>Halomicrobium</i> sp. ZPS1 — показательный случай. В этом геноме два гена 16S рРНК, поэтому энтропия равна количеству нуклеотидов, которые не совпадают в их последовательностях. Вычислив процент идентичности (сокращённо ПИ) последовательностей этих двух генов, получим 90,70%! Напомним, что за порог ПИ последовательностей генов 16S рРНК для отнесения организмов да отдельных <em>родов</em> обычно принимают (впрочем, произвольно) значение 95%.</p>
<hr>

<div id="primers-coverages" class="pad-anchor"></div>
<h3>Спектр действия<sup>*</sup> пар праймеров к разным V-участкам бактериальных генов 16S рРНК</h3>

<p><sup>*</sup> Спектр действия пары праймеров — это процент геномов, которые содержат как минимум один ген 16S рРНК, который возможно амплифицировать с помощью ПЦР используя эту пару праймеров. Подробности описаны в нашей <a href="https://doi.org/10.1016/j.resmic.2022.103936">статье про RiboGrove</a>.</p>
<p>В таблицах ниже можно найти спектры действия пар праймеров, которые обычно используются для амплификации бактериальных и архейных генов 16S рРНК („бактериальные“ и „архейные“ праймеры).</p>
<p>Есть более подробная таблица — в файле <span class="samp">primer_pair_genomic_coverage.tsv</span> в метаданных <a href="#downloads">metadata</a>. Та таблица содержит спектры действия не только для отделов (phyla), но и для каждого класса, порядка, семейства, рода и вида. Кроме того, та таблица содержит спектры действия для пары праймеров 1115F–1492R (участок V7–V9). В таблицах ниже нет этой пары, чтобы они не были слишком громоздкими.</p>

<table class="sum-table"><caption>Бактерияльные гены, „бактерияльные“ праймеры</caption>
<tbody class="primer-cov-tbody">
<tr>
  <th class="alnleft" rowspan="2">Отдел</th>
  <th class="numcol" rowspan="2">Количество<br />геномов</th>
  {% for v_region_name in bacterial_primer_pairs.values() %}
  <th class="numcol">{{ v_region_name.replace('-', '–') }}</th>
  {% endfor %}
</tr>
<tr>
  {% for primer_pair_name in bacterial_primer_pairs.keys() %}
  <th class="numcol">{{ primer_pair_name.replace('-', '–') }}<br />(%)</th>
  {% endfor %}
</tr>
{% for _, row in ribogrove_primers_cov_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
  <td>{{ italicize_candidatus(row['Phylum']) }}</td>
  <td class="numcol">{{ row['num_genomes'] }}</td>
  {% for primer_pair_name in bacterial_primer_pairs.keys() %}
  <td class="numcol">{{ row[primer_pair_name] }}</td>
  {% endfor %}
</tr>
{% endfor %}
</tbody>
</table>
<br />

<table class="sum-table"><caption>Архейные гены, „архейные“ праймеры</caption>
<tbody class="primer-cov-tbody">
<tr>
  <th class="alnleft" rowspan="2">Отдел</th>
  <th class="numcol" rowspan="2">Количество<br />геномов</th>
  {% for v_region_name in archaeal_primer_pairs.values() %}
  <th class="numcol">{{ v_region_name.replace('-', '–') }}</th>
  {% endfor %}
</tr>
<tr>
  {% for primer_pair_name in archaeal_primer_pairs.keys() %}
  <th class="numcol">{{ primer_pair_name.replace('-', '–') }}<br />(%)</th>
  {% endfor %}
</tr>
{% for _, row in ribogrove_primers_cov_df.query('Domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
  <td>{{ italicize_candidatus(row['Phylum']) }}</td>
  <td class="numcol">{{ row['num_genomes'] }}</td>
  {% for primer_pair_name in archaeal_primer_pairs.keys() %}
  <td class="numcol">{{ row[primer_pair_name] }}</td>
  {% endfor %}
</tr>
{% endfor %}
</tbody>
</table>
<br />

<details>
<summary><b>Бактериальные гены, „архейные“ праймеры</b></summary>
<div id="primers-coverages" class="pad-anchor"></div>
<table class="sum-table"><caption>Бактериальные гены, „архейные“ праймеры</caption>
<tbody class="primer-cov-tbody">
<tr>
  <th class="alnleft" rowspan="2">Отдел</th>
  <th class="numcol" rowspan="2">Количество<br />геномов</th>
  {% for v_region_name in archaeal_primer_pairs.values() %}
  <th class="numcol">{{ v_region_name.replace('-', '–') }}</th>
  {% endfor %}
</tr>
<tr>
  {% for primer_pair_name in archaeal_primer_pairs.keys() %}
  <th class="numcol">{{ primer_pair_name.replace('-', '–') }}<br />(%)</th>
  {% endfor %}
</tr>
{% for _, row in ribogrove_primers_cov_df.query('Domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
  <td>{{ italicize_candidatus(row['Phylum']) }}</td>
  <td class="numcol">{{ row['num_genomes'] }}</td>
  {% for primer_pair_name in archaeal_primer_pairs.keys() %}
  <td class="numcol">{{ row[primer_pair_name] }}</td>
  {% endfor %}
</tr>
{% endfor %}
</tbody>
</table>
</details>
<br />

<details>
<summary><b>Архейные гены, „бактериальные“ праймеры</b></summary>
<table class="sum-table"><caption>Архейные гены, „бактериальные“ праймеры</caption>
  <tbody class="primer-cov-tbody">
  <tr>
    <th class="alnleft" rowspan="2">Отдел</th>
    <th class="numcol" rowspan="2">Количество<br />геномов</th>
    {% for v_region_name in bacterial_primer_pairs.values() %}
    <th class="numcol">{{ v_region_name.replace('-', '–') }}</th>
    {% endfor %}
  </tr>
  <tr>
    {% for primer_pair_name in bacterial_primer_pairs.keys() %}
    <th class="numcol">{{ primer_pair_name.replace('-', '–') }}<br />(%)</th>
    {% endfor %}
  </tr>
  {% for _, row in ribogrove_primers_cov_df.query('Domain == "Archaea"').iterrows() %}
  <tr class="sumtab-row">
    <td>{{ italicize_candidatus(row['Phylum']) }}</td>
    <td class="numcol">{{ row['num_genomes'] }}</td>
    {% for primer_pair_name in bacterial_primer_pairs.keys() %}
    <td class="numcol">{{ row[primer_pair_name] }}</td>
    {% endfor %}
  </tr>
  {% endfor %}
  </tbody>
</table>
</details>
<br />


<details>
<summary><b>Праймеры, спектр действия которых оценивался</b></summary>
<table class="sum-table"><caption>Праймеры, спектр действия которых оценивался</caption>
<tbody>
  <tr><th>Название праймера</th><th>Последовательность</th><th>Ссылка</th></tr>
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
  <tr class="sumtab-row"><td>SSU1ArF</td><td>TCCGGTTGATCCYGCBRG</td><td><a href="https://doi.org/10.1111/1758-2229.12684">Bahram et al., 2018</a></td></tr>
  <tr class="sumtab-row"><td>SSU520R</td><td>GCTACGRRYGYTTTARRC</td><td><a href="https://doi.org/10.1111/1758-2229.12684">Bahram et al., 2018</a></td></tr>
  <tr class="sumtab-row"><td>340f</td><td>CCCTAYGGGGYGCASCAG</td><td><a href="https://doi.org/10.1016/j.mimet.2010.10.001">Gantner et al., 2011</a></td></tr>
  <tr class="sumtab-row"><td>806rB</td><td>GGACTACNVGGGTWTCTAAT</td><td><a href="https://doi.org/10.3354/ame01753 ">Appril et al., 2015</a></td></tr>
  <tr class="sumtab-row"><td>349f</td><td>GYGCASCAGKCGMGAAW</td><td><a href="https://doi.org/10.1128/AEM.66.11.5066-5072.2000">Takai and Horikoshi, 2000</a></td></tr>
  <tr class="sumtab-row"><td>519r</td><td>TTACCGCGGCKGCTG</td><td><a href="https://doi.org/10.1093/nar/gks808">Klindworth et al., 2013</a></td></tr>
  <tr class="sumtab-row"><td>515fB</td><td>GTGYCAGCMGCCGCGGTAA</td><td><a href="https://doi.org/10.1111/1462-2920.13023">Parada et al., 2015</a></td></tr>
  <tr class="sumtab-row"><td>Parch519f</td><td>CAGCCGCCGCGGTAA</td><td><a href="https://doi.org/10.1128/aem.63.9.3367-3373.1997">Ovreås et al., 1997</a></td></tr>
  <tr class="sumtab-row"><td>Arch915r</td><td>GTGCTCCCCCGCCAATTCCT</td><td><a href="https://doi.org/10.1128/aem.60.4.1232-1240.1994">Raskin et al., 1994</a></td></tr>
  <tr class="sumtab-row"><td>1106F</td><td>TTWAGTCAGGCAACGAGC</td><td><a href="https://doi.org/10.1016/j.soilbio.2007.05.030">Watanabe et al., 2007</a></td></tr>
  <tr class="sumtab-row"><td>Ar1378R<sup> **</sup></td><td>TGTGCAAGGAGCAGGGAC</td><td><a href="https://doi.org/10.1016/j.soilbio.2007.05.030">Watanabe et al., 2007</a></td></tr>
  <tr class="sumtab-row"><td>A751F</td><td>CCGACGGTGAGRGRYGAA</td><td><a href="https://doi.org/10.1016/j.mimet.2003.08.009">Baker et al., 2003</a></td></tr>
  <tr class="sumtab-row"><td>SSU1492Rngs</td><td>CGGNTACCTTGTKACGAC</td><td><a href="https://doi.org/10.1111/1758-2229.12684">Bahram et al., 2018</a></td></tr>
  <tr class="sumtab-row"><td>SSU280ArR</td><td>TCAGWNYCCNWCTCSRGG</td><td><a href="https://doi.org/10.1111/1758-2229.12684">Bahram et al., 2018</a></td></tr>
  <tr class="sumtab-row"><td>SSU470R</td><td>DCNGCNGGTDTTACCGCG</td><td><a href="https://doi.org/10.1111/1758-2229.12684">Bahram et al., 2018</a></td></tr>
  <tr class="sumtab-row"><td>SSU468R</td><td>GNDCNGCNGGTDTTACCG</td><td><a href="https://doi.org/10.1111/1758-2229.12684">Bahram et al., 2018</a></td></tr>
  <tr class="sumtab-row"><td>A519R</td><td>GGTDTTACCGCGGCKGCTG</td><td><a href="https://doi.org/10.1371/journal.pone.0007401">Wang and Qian, 2009</a></td></tr>
  <tr class="sumtab-row"><td>SSU666ArR</td><td>HGCYTTCGCCACHGGTRG</td><td><a href="https://doi.org/10.1111/1758-2229.12684">Bahram et al., 2018</a></td></tr>
  <tr class="sumtab-row"><td>SSU1000ArR</td><td>GGCCATGCAMYWCCTCTC</td><td><a href="https://doi.org/10.1111/1758-2229.12684">Bahram et al., 2018</a></td></tr>
  <tr class="sumtab-row"><td>UA1204R</td><td>TTMGGGGCATRCIKACCT</td><td><a href="https://doi.org/10.1016/j.mimet.2003.08.009">Baker et al., 2003</a></td></tr>
</tbody>
</table>
<p><sup>*</sup> Праймеры 341F и 785R используются в <a href="https://support.illumina.com/downloads/16s_metagenomic_sequencing_library_preparation.html">протоколе</a> подготовки библиотеки для секвенирования ампликонов участков V3–V4 генов 16S рРНК на приборе Illumina MiSeq.</p>
<p><sup>**</sup> Исходное название праймера Ar1378R — 1378R. Мы изменили название, чтобы избежать путаницы.</p>
</details>


<div id="searching-data" class="pad-anchor"></div>
<h2>Поиск данных в базе данных RiboGrove</h2>
<p>RiboGrove — очень минималистичная база данных: она представляет собой несколько файлов в формате fasta, а также метаданные. Поэтому, расширенные средства поиска для неё недоступны. Мы признаём эту проблему и оставляем рекомендации по поиску данных в RiboGrove. Рекомендации помогут вам изучать и делать выборки данных из БД RiboGrove.</p>
<div id="header-format" class="pad-anchor"></div>
<h3>Формат fasta-заголовка</h3>
<p>Fasta-файлы, составляющие базу данных RiboGrove имеют следующий формат заголовка:</p>
<p class="samp-highl samp-vwide" style="font-size:14px">&gt;GCF_000978375.1:NZ_CP009686.1:8908-10459:plus ;d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__cereus; category:1</p>
<p>Основные части заголовка отделены друг от друга пробелами. Заголовок состоит из трёх таких частей:</p>
<ol>
<li>Идентификатор последовательности (seqID): <span class="samp-highl">GCF_000978375.1:NZ_CP009686.1:8908-10459:plus</span>. SeqID, в свою очередь, состоит из четырёх частей, отделённых друг от друга двоеточиями (<span class="samp-highl">:</span>):
<ol>
<li>Код доступа сборки генома, которому принадлежит последовательность гена: <span class="samp-highl">GCF_000978375.1</span>.</li>
<li>Код доступа геномной последовательности из БД RefSeq, которой принадлежит последовательность гена: <span class="samp-highl">NZ_CP009686.1</span>.</li>
<li>Координаты гена в пределах геномной последовательности из RefSeq: <span class="samp-highl">8908-10459</span> (отсчёт координат начинается с единицы, правая и левая координаты включены в интервал).</li>
<li>Знак цепи геномной последовательности из RefSeq, на которой расположен ген: <span class="samp-highl">plus</span> (либо <span class="samp-highl">minus</span>).</li>
</ol>
</li>
<li>Строка таксономии, а именно названия домена (<span class="samp-highl">Bacteria</span>), отдела (<span class="samp-highl">Firmicutes</span>), класса (<span class="samp-highl">Bacilli</span>), порядка (<span class="samp-highl">Bacillales</span>), семейства (<span class="samp-highl">Bacillaceae</span>) и рода (<span class="samp-highl">Bacillus</span>), а также видовой эпитет (<span class="samp-highl">cereus</span>).<br>
Перед каждым названием пишется префикс, который обозначает ранг: <span class="samp-highl">d__</span> для домена, <span class="samp-highl">p__</span> для отдела, <span class="samp-highl">c__</span> для класса, <span class="samp-highl">o__</span> для порядка, <span class="samp-highl">f__</span> для семейства, <span class="samp-highl">g__</span> для рода и <span class="samp-highl">s__</span> для видового эпитета. В префиксах после буквы пишутся <strong>два</strong> знака нижнего подчёркивания.<br>
Названия разделены и окружены точками с запятой (<span class="samp-highl">;</span>).</li>
<li>Категория генома, и которого происходит ген: (<span class="samp-highl">category:1</span>).</li>
</ol>
<h3>Выборка последовательностей</h3>
<p>Отбирать последовательности из fasta-файлов можно с помощью программы <strong>Seqkit</strong> (<a href="https://github.com/shenwei356/seqkit">ГитХаб-репозиторий</a>, <a href="https://bioinf.shenwei.me/seqkit/">документация</a>). Это бесплатная, кроссплатформенная, многофункциональная и весьма шустрая программа; она может обрабатывать gzipʼнутые и простые fasta-файлы. Для выборки данных из fasta-файлов удобно использовать программы <span class="samp">seqkit grep</span> и <span class="samp">seqkit seq</span>.</p>
<div id="subset-by-header" class="pad-anchor"></div>
<h4>Выборка последовательностей по заголовку</h4>
<p>Пусть мы делаем выборку из файла <span class="samp">ribogrove_6.212_sequences.fasta.gz</span>. Полезными могут быть следующие примеры команд для выборки последовательностей с помощью <span class="samp">seqkit grep</span>:</p>
<p><strong>Пример 1</strong>. Выбрать индивидуальную последовательность по идентификатору (seqID).</p>
<p class="samp-highl samp-vwide">seqkit grep -p "GCF_000978375.1:NZ_CP009686.1:8908-10459:plus" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опция <span class="samp">-p</span> задаёт строку, которую программа будет искать в fasta-заголовках (на самом дела, только в той части заголовков, в которой записан seqID).</p>
<p><strong>Пример 2</strong>. Выбрать все последовательности генов, которые происходят из индивидуальной геномной последовательности из БД RefSeq, по коду доступа геномной последовательности <span class="samp">NZ_CP009686.1</span>.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ":NZ_CP009686.1:" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Тут, нужно указывать ещё две опции: <span class="samp">-n</span> и <span class="samp">-r</span>. Первая указывает программе, что ей следует искать вхождение строки-запроса в заголовках целиком, а не только в seqID. Вторая опция указывает программе, что ей следует искать не только полные совпадения (строка-запрос идентична полному заголовку), но и частичные (когда строка-запрос является подстрокой заголовка).</p>
<p class="samp-comment">Чтобы обеспечить специфичность поиска, окружайте код доступа двоеточиями (<span class="samp-highl">:</span>).</p>
<p><strong>Пример 3</strong>. Выбрать все последовательности генов индивидуального генома (код доступа геномной сборки <a href="https://www.ncbi.nlm.nih.gov/assembly/GCF_019357495.1">GCF_019357495.1</a>).</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "GCF_019357495.1:" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Чтобы обеспечить специфичность поиска, пишите двоеточие (<span class="samp-highl">:</span>) после кода доступа сборки.</p>
<p><strong>Пример 4</strong>. Выбрать все последовательности генов актинобактерий.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ";p__Actinobacteria;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Чтобы обеспечить специфичность поиска, стоит окружать названия таксонов точками с запятой (<span class="samp-highl">;</span>).</p>
<p><strong>Пример 5</strong>. Выбрать все последовательности генов, которые происходят из геномов 1-й категории.</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "category:1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Пример 6</strong>. Выбрать все последовательности генов, кроме тех, которые принадлежат фирмикутам.</p>
<p class="samp-highl samp-vwide">seqkit grep -nvrp ";p__Firmicutes;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Обратите внимание на опцию <span class="samp">-v</span> внутри строки с опциями (<span class="samp">-nvrp</span>). Эта опция инвертирует поведение программы, и теперь она находит только те последовательности, в заголовках которых нет подстроки «<span class="samp">;p__Firmicutes;</span>».</p>
<div id="subset-by-len" class="pad-anchor"></div>
<h4>Выборка последовательностей по длине</h4>
<p>Отбирать последовательности по длине можно с помощью программы <span class="samp">seqkit seq</span>.</p>
<p><strong>Пример 1</strong>. Выбрать все последовательности длиннее 1600 п.н.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1601 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опция <span class="samp">-m</span> указывает программе минимальную длину последовательности, которую программа подаст на выход.</p>
<p><strong>Пример 2</strong>. Выбрать все последовательности короче 1500 п.н.</p>
<p class="samp-highl samp-vwide">seqkit seq -M 1499 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опция <span class="samp">-M</span> указывает программе максимальную длину последовательности, которую программа подаст на выход.</p>
<p><strong>Пример 3</strong>. Выбрать все последовательности, длина которых находится в промежутке [1500, 1600] п.н.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1500 -M 1600 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<div id="select-head" class="pad-anchor"></div>
<h4>Извлечение информации из заголовков последовательностей</h4>
<p>Это можно сделать с помощью программы <span class="samp">seqkit seq</span>.</p>
<p><strong>Пример 1</strong>. Выбрать все заголовки.</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опция <span class="samp">-n</span> указывает программе подавать на выход только fasta-заголовки.</p>
<p><strong>Пример 2</strong>. Выбрать все идентификаторы последовательностей (части заголовков до первого пробела).</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опция <span class="samp">-i</span> указывает программе подавать на выход только идентификаторы последовательностей.</p>
<p><strong>Пример 3</strong>. Выбрать все коды доступа геномных последовательностей из базы данных RefSeq.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f2 -d':' | sort | uniq</p>
<p class="samp-comment">Такая команда сработает успешно, если в вашей системе установлены программы <span class="samp">cut</span>, <span class="samp">sort</span> и <span class="samp">uniq</span> (в системах Linux и Mac OS эти программы обычно предустановлены).</p>
<p><strong>Пример 4</strong>. Выбрать все коды доступа геномных сборок.</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | cut -f1 -d':' | sort | uniq</p>
<p class="samp-comment">Такая команда сработает успешно, если в вашей системе установлены программы <span class="samp">cut</span>, <span class="samp">sed</span>, <span class="samp">sort</span> и <span class="samp">uniq</span> (в системах Linux и Mac OS эти программы обычно предустановлены).</p>
<p><strong>Пример 5</strong>. Выбрать все названия отделов организмов.</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | grep -Eo ';p__[^;]+' | sed -E 's/;|p__//g' | sort | uniq</p>
<p class="samp-comment">Такая команда сработает успешно, если в вашей системе установлены программы <span class="samp">grep</span>, <span class="samp">sed</span>, <span class="samp">sort</span> и <span class="samp">uniq</span> (в системах Linux и Mac OS эти программы обычно предустановлены).</p>
{% if not archive %}
<hr>
<div id="contacts" class="pad-anchor"></div>
<h2>Контакты</h2>
<p>Если у вас есть вопрос про RiboGrove, пожалуйста, обращайтесь к Максиму Сиколенко на эл. адрес sikolenko<img class="sabaka" src="/wp-content/uploads/cager/pes2.gif" alt="[ at ]" align="bottom">mbio.bas-net.by или maximdeynonih<img class="sabaka" src="img/pes2.gif" alt="[ at ]" align="bottom">gmail.com.</p>
<hr>
<div id="citing-ribogrove" class="pad-anchor"></div>
<h2>Цитирование RiboGrove</h2>
<p>Если RiboGrove была полезна вам в вашей работе, пожалуйста, цитируйте следующую публикацию:</p>
<p>Maxim A. Sikolenko, Leonid N. Valentovich. “RiboGrove: a database of full-length prokaryotic 16S rRNA genes derived from completely assembled genomes” // Research in Microbiology, Volume 173, Issue 4, May 2022, 103936.<br> (ЦИО: <a href="https://doi.org/10.1016/j.resmic.2022.103936">10.1016/j.resmic.2022.103936</a>).</p>
<hr>
<div id="faq" class="pad-anchor"></div>
<h2>Вопросы, которые возникают у людей про RiboGrove</h2>
<div><details>
  <summary><b>1. Как создать файл с таксономией, с которым смогли бы работать программы из пакета QIIME2?</b></summary>
  <p>Для этого вы можете воспользоваться скриптом <a href="https://github.com/masikol/ribogrove-tools/blob/main/data_transformation_scripts/make_qiime_taxonomy_file.py"><span class="samp">make_qiime_taxonomy_file.py</span></a>: он трансформирует файл <span class="samp">metadata/taxonomy.tsv</span> в файл, с которым смогут работать программы из QIIME2. Справку-руководство для этого скрипта можете найти в соответствующем <a href="https://github.com/masikol/ribogrove-tools/blob/main/data_transformation_scripts/README.md">файле README</a>.</p>
</details></div><br>
<div><details>
  <summary><b>2. Как сохранить в файл данные, отобранные с помощью Seqkit?</b></summary>
  <p>В соответствующей дискуссии на форуме уже дали несколько полезных ответов и рекомендаций: <a href="https://bioinformatics.stackexchange.com/questions/20915/how-do-i-save-selected-sequences-in-seqkit-to-a-file">https://bioinformatics.stackexchange.com/questions/20915/how-do-i-save-selected-sequences-in-seqkit-to-a-file</a>.</p>
</details></div><br>
<div><details>
  <summary><b>3. Как найти нужные записи в fasta-файле по последовательности с помощью Seqkit?</b></summary>
  <p>В соответствующей дискуссии на форуме уже дали несколько полезных ответов и рекомендаций: <a href="https://www.biostars.org/p/9561418">https://www.biostars.org/p/9561418</a>.</p>
</details></div>
{% endif %}