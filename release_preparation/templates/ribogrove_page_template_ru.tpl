<img src="/cager/img/RiboGrove_logo.gif"><br /><br />
<h2>Содержание</h2>
<ul>
<li><a href="#overview">Описание</a></li>
<li><a href="#downloads">Файлы</a></li>
<li><a href="#summary">Сводка</a>
<ul style="list-style-type: none;">
<li><a href="#db-size">Размер базы данных RiboGrove</a></li>
<li><a href="#gene-lens">Длина генов 16S рРНК</a></li>
<li><a href="#copy-number">Число копий генов 16S рРНК</a></li>
<li><a href="#top-longest">Топ-10 самых длинных генов 16S рРНК</a></li>
<li><a href="#top-shortest">Топ-10 самых коротких генов 16S рРНК</a></li>
<li><a href="#top-copy-num">Tоп-10 геномов с наибольшим числом копий генов 16S рРНК</a></li>
<li><a href="#top-var">Tоп-10 геномов с наибольшей внутригеномной изменчивостью генов 16S рРНК</a></li>
</ul>
</li>
<li><a href="#searching-data">Поиск данных в базе данных RiboGrove</a>
<ul style="list-style-type: none;">
<li><a href="#subset-by-header">Выборка последовательностей по заголовку</a></li>
<li><a href="#subset-by-len">Выборка последовательностей по длине</a></li>
<li><a href="#select-head">Извлечение информации из заголовков последовательностей</a></li>
</ul>
</li>
</ul>
<hr />
<h2>Описание</h2>
<div id="overview" class="pad-anchor"></div>
<p>RiboGrove – это база данных последовательностей генов 16S рРНК бактерий и архей. RiboGrove содержит исключительно полноразмерные последовательнсти генов 16S рРНК, т.к. при составлении базы данных не использовались данные секвенирования ампликоннов. Другими, словами, в отличие от сходных информационных ресурсов, посвящённых последовательностям (генов) рРНК – таких как <a href="https://www.arb-silva.de/">Silva</a> или <a href="http://rdp.cme.msu.edu/">RDP</a> – RiboGrove не содержит неполных последовательностей генов (например, только участков V3-V4).</p>
<p>Для составления базы данных RiboGrove использовались полностью собранные геномы прокариот, депонированные в базу данных <a href="https://www.ncbi.nlm.nih.gov/refseq/">RefSeq</a>.</p>
<p>Программы, с помощью которых RiboGrove была создана, доступны в репозитории на ГитХабе: <a href="https://github.com/masikol/ribogrove-paper">ribogrove-paper</a>.</p>
<hr />
<h2>Файлы</h2>
<div id="downloads" class="pad-anchor"></div>
<h3>Актуальная версия &ndash; {{ ribogrove_release_number }} ({{ ribogrove_release_date }})</h3>
<p>Выпуск RiboGrove {{ ribogrove_release_number }} основан на данных, полученных из базы данных RefSeq версии 207.</p>
<ul>
<li>Fasta-файл полноразмерных последовательностей генов 16S рРНК. <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Скачать (gzip'нутый fasta-файл, {{ final_fasta_fsize_fmt }} Мб)</a></li>
<li>"Сырая" версия fasta-файла, указанного выше. Данный файл содержит в том числе и неполные последовательности генов. <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/raw_ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz">Скачать (gzip'нутый fasta-файл, {{ raw_fasta_fsize_fmt }} Мб)</a></li>
<li>Метаданные <a href="/cager/ribogrove_releases/{{ ribogrove_release_number }}/metadata_ribogrove_{{ ribogrove_release_number }}.zip">Скачать (zip-архив {{ metadata_fsize_fmt }} Мб)</a></li>
</ul>
<p>Fasta-файлы сжаты с помощью программы gzip, а метаданные – с помощью программы zip. Чтобы разархивировать их, пользователи Linux и Mac OS могут воспользоваться программами gzip и zip, соответственно; эти программы должны быть предустановлены. Пользователи же Windows могут воспользоваться бесплатной программой <a href="https://www.7-zip.org/">7-Zip</a>.</p>
<h3>Архив выпусков RiboGrove</h3>
<p>Вы можете найти все выпуски RiboGrove в <a href="/cager/ru/content/60-ribogrove-archive">архиве выпусков RiboGrove</a>.</p>
<hr />
<h2>Сводка</h2>
<div id="summary" class="pad-anchor"></div>
<table class="sum-table"><caption>Размер базы данных RiboGrove
<div id="db-size" class="pad-anchor"></div>
</caption>
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
<table class="sum-table"><caption>Длина генов 16S рРНК
<div id="gene-lens" class="pad-anchor"></div>
</caption>
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
<p><sup>*</sup> &ndash; метрики, помеченные звёздочкой, были расчитаны с предварительной нормализацией данных, а именно: медианные значения длины внутри каждого вида использовались для расчёта этих метрик.</p>
<table class="sum-table"><caption>Число копий генов 16S рРНК
<div id="copy-number" class="pad-anchor"></div>
</caption>
<tbody>
<tr><th class="numcol">Число копий <sup>*</sup></th><th class="numcol">Количество видов</th><th class="numcol">Процент видов, %</th></tr>
{% for _, row in ribogrove_copy_number_df.iterrows() %}
<tr class="sumtab-row">
<td class="numcol">{{ row['copy_number'] }}</td>
<td class="numcol">{{ row['number_of_species'] }}</td>
<td class="numcol">{{ row['percent_of_species'] }}</td>
</tr>
{% endfor %}
</tbody>
</table>
<p><sup>*</sup> &ndash; медианное число копий в пределах вида.</p>

<table class="sum-table"><caption>Топ-10 самых длинных генов 16S рРНК
<div id="top-longest" class="pad-anchor"></div>
</caption>
<tbody>
<tr><th class="alnleft">ID геномной сборки</th><th class="numcol">Длина гена, п.н.</th><th class="alnleft">Идентификатор(ы) послед-ти(ей)<br />в БД RiboGrove</th><th class="alnleft">Организм</th></tr>
<tr>
<td colspan="4" class="subhead">Бактрии</td>
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
<td colspan="4" class="subhead">Археи</td>
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

<table class="sum-table"><caption>Топ-10 самых коротких генов 16S рРНК
<div id="top-shortest" class="pad-anchor"></div>
</caption>
<tbody>
<tr><th class="alnleft">ID геномной сборки</th><th class="numcol">Длина гена, п.н.</th><th class="alnleft">Идентификатор(ы) послед-ти(ей)<br />в БД RiboGrove</th><th class="alnleft">Организм</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерии</td>
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
<td colspan="4" class="subhead">Археи</td>
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

<table class="sum-table"><caption>Tоп-10 геномов с наибольшим числом копий генов 16S рРНК
<div id="top-copy-num" class="pad-anchor"></div>
</caption>
<tbody>
<tr><th class="alnleft">ID геномной сборки</th><th class="numcol">Число копий</th><th class="alnleft">Организм</th></tr>
<tr>
<td colspan="4" class="subhead">Бактерии</td>
</tr>
{% for _, row in ribogrove_top_copy_numbers_df.query('domain == "Bacteria"').iterrows() %}
    <tr class="sumtab-row">
    <td><a href="https://ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
    <td class="numcol">{{ row['copy_number'] }}</td>
    <td><i>{{ row['species'] }}</i></td>
    </tr>
{% endfor %}
<tr>
<td colspan="4" class="subhead">Археи</td>
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

<table class="sum-table"><caption>Tоп-10 геномов с наибольшей внутригеномной изменчивостью генов 16S рРНК
<div id="top-var" class="pad-anchor"></div>
</caption>
<tbody>
<tr><th class="alnleft">ID геномной сборки</th><th class="numcol">Сумма энтропии<sup>*</sup>, биты</th><th class="numcol">Средняя энтропия<sup>*</sup>, биты</th><th class="numcol">Количество изменсивых позиций</th><th class="numcol">Число копий</th><th class="alnleft">Организм</th></tr>

<tr>
<td colspan="6" class="subhead">Бактерии</td>
</tr>
{% for _, row in ribogrove_top_intragenomic_var_df.query('domain == "Bacteria"').iterrows() %}
<tr class="sumtab-row">
<td><a href="https://www.ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
<td class="numcol">{{ row['sum_entropy'] }}</td>
<td class="numcol">{{ row['mean_entropy'] }}</td>
<td class="numcol">{{ row['num_var_cols'] }}</td>
<td class="numcol">{{ row['copy_numer'] }}</td>
<td><i>{{ row['species'] }}</i></td>
</tr>
{% endfor %}
<tr>
<td colspan="6" class="subhead">Археи</td>
</tr>
{% for _, row in ribogrove_top_intragenomic_var_df.query('domain == "Archaea"').iterrows() %}
<tr class="sumtab-row">
<td><a href="https://www.ncbi.nlm.nih.gov/assembly/{{ row['ass_id'] }}">{{ row['ass_id'] }}</a></td>
<td class="numcol">{{ row['sum_entropy'] }}</td>
<td class="numcol">{{ row['mean_entropy'] }}</td>
<td class="numcol">{{ row['num_var_cols'] }}</td>
<td class="numcol">{{ row['copy_numer'] }}</td>
<td><i>{{ row['species'] }}</i></td>
</tr>
{% endfor %}
</tbody>
</table>

<p><sup>*</sup> &ndash; энтропия в данном случае это <a href="https://ru.wikipedia.org/wiki/Информационная_энтропия">энтропия Шенона</a>, вычисленная для каждого столбца множественного выравнивания полноразмерных последовательностей генов 16S рРНК индивидуального генома. Вычисленная энтропия далее суммировалась (столбец таблицы "Сумма энтропии") или усреднялась (столбец таблицы "Средняя энтропия").</p>
<hr />
<h2>Поиск данных в базе данных RiboGrove</h2>
<div id="searching-data" class="pad-anchor"></div>
<p>RiboGrove &ndash; очень минималистичная база данных: она представляет собой несколько файлов в формате fasta, а также метаданные. Поэтому, расширенные средства поиска для неё недоступны. Мы признаём эту проблему и оставляем рекоммендации по поиску данных в RiboGrove. Рекоммендации помогут вам изучать и делать выборки данных из БД RiboGrove.</p>
<h3>Формат fasta-заголовка</h3>
<p>Fasta-файлы, составляющие базу данных RiboGrove имеют следующий формат заголовка:</p>
<p class="samp-highl samp-vwide">&gt;NZ_CP079719.1:86193-87742_plus Bacillus_velezensis ;Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus; category:2</p>
<p>Основные части заголовка отделены друг от друга пробелами. Заголовок состоит из четырёх таких частей:</p>
<ol>
<li>Идентификатор последовательности (seqID): <span class="samp-highl">NZ_CP079719.1:86193-87742_plus</span>. SeqID, в свою очередь, состоит из трёх частей:<ol>
<li>код доступа геномной последовательности из БД RefSeq, которой принадлежит последовательность гена: <span class="samp-highl">NZ_CP079719.1</span>;</li>
<li>координаты гена в пределах геномной последовательности из RefSeq: <span class="samp-highl">86193-87742</span> (отсчёт координат начинается с единицы, правая и левая координаты включительны);</li>
<li>знак цепи геномной последовательности из RefSeq, на которой расположен ген: <span class="samp-highl">plus</span>.</li>
</ol></li>
<li>Название организма, в соответствии с базой данных <a href="https://www.ncbi.nlm.nih.gov/taxonomy">NCBI Taxonomy</a>: <span class="samp-highl">Bacillus_velezensis</span>.</li>
<li>Строка таксономии, а именно названия домена (<span class="samp-highl">Bacteria</span>), отдела (<span class="samp-highl">Firmicutes</span>), класса (<span class="samp-highl">Bacilli</span>), порядка (<span class="samp-highl">Bacillales</span>), семейства (<span class="samp-highl">Bacillaceae</span>), и рода (<span class="samp-highl">Bacillus</span>). Названия разделены и окружены точками с запятой (<span class="samp-highl">;</span>).</li>
<li>Категория генома, и которого происходят ген: (<span class="samp-highl">category:2</span>. Сборки геномов 1-й категории наиболее достоверны, а геномов 3-й категории   наименее достоверны.</li>
</ol>
<h3>Выборка последовательностей</h3>
<p>Вы можете делать выборки последовательностей из fasta-файлов с помощью программы <strong>Seqkit</strong> (<a href="https://github.com/shenwei356/seqkit">ГитХаб-репозиторий</a>, <a href="https://bioinf.shenwei.me/seqkit/">документация</a>). Это бесплатная, кроссплатформенная, многофункциональная и весьма шустрая программа; она может обрабатывать gzip'нутые и простые fasta-файлы. А конкретно, для выборки данных из fasta-файлов следует использовать (под)программу <span class="samp">seqkit grep</span>  .</p>
<h4>Выборка последовательностей по заголовку</h4>
<div id="subset-by-header" class="pad-anchor"></div>
<p>Пусть мы делаем выборку из файла <span class="samp">ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</span>. Следующие примеры команд для выборки последовательностей с помощью <span class="samp">seqkit grep</span> могут быть вам полезны:</p>
<p><strong>Пример 1</strong>. Выбрать индивидуальную последовательность по идентификатору (seqID)</p>
<p class="samp-highl samp-vwide">seqkit grep -p "NZ_CP079719.1:86193-87742_plus" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опция <span>-p</span> задаёт строку, которую программа будет искать в fasta-заголовках (на самом дела, только в той части заголовков, в которой записан seqID).</p>
<p><strong>Пример 2</strong>. Выбрать все последовательности генов, которые происходят из индивидуальной геномной последовательности single из БД RefSeq, по коду доступа геномной последовательности <span class="samp">NZ_CP079719.1</span></p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "NZ_CP079719.1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Тут, нужно указывать ещё две опции: <span class="samp">-n</span> и <span class="samp">-r</span>. Первая указывает программе, что ей следует искать вхождение строки-запроса в заголовках целиком, а не только в seqID. Вторая опция указывает программе, что ей следует искать не только полные совпадения (строка-запрос идентична полному заголовку), но и частичные (когда строка-запрос является подстрокой заголовка).</p>
<p><strong>Пример 3</strong>. Выбрать все последовательности генов индивидуального генома (ID геномной сборки <a href="https://www.ncbi.nlm.nih.gov/assembly/10577151">10577151</a>), который состоит из двух репликонов: <a href="https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP079110.1">NZ_CP079110.1</a> и <a href="https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP079111.1">NZ_CP079111.1</a></p>
<p class="samp-highl samp-vwide">seqkit grep -nr -p "NZ_CP079110.1" -p "NZ_CP079111.1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Пример 4</strong>. Выбрать все последовательности генов актинобактерий</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp ";Actinobacteria;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">На всякий случай стоит окружать названия таксонов точками с запятой (<span class="samp-highl">;</span>).</p>
<p><strong>Пример 5</strong>. Выбрать все последовательности генов, которые происходят из геномой 1-й категории</p>
<p class="samp-highl samp-vwide">seqkit grep -nrp "category:1" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p><strong>Пример 6</strong>. Выбрать все последовательности генов, кроме тех, которые принадлежат фирмикутам</p>
<p class="samp-highl samp-vwide">seqkit grep -nvrp ";Firmicutes;" ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Обратите внимание на опцию <span class="samp">-v</span> внутри строки с опциями (<span class="samp">-nvrp</span>). Эта опция инвертирует поведение программы, и теперь она находит только те последовательности, в заголовках которых нет подстроки <span class="samp">;Firmicutes;</span>.</p>
<h4>Выборка последовательностей по длине</h4>
<div id="subset-by-len" class="pad-anchor"></div>
<p>Пусть мы делаем выборку из файла <span class="samp">ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</span>. Следующие примеры команд для выборки последовательностей по длине с помощью <span class="samp">seqkit seq</span> могут быть вам полезны:</p>
<p><strong>Пример 1</strong>. Выбрать все последовательности длиннее 1′600 п.н.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1601 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опция <span class="samp">-m</span> указывает программе минимальную длину последовательности, которую программа подаст на выход.</p>
<p><strong>Пример 2</strong>. Выбрать все последовательности короче 1′500 п.н.</p>
<p class="samp-highl samp-vwide">seqkit seq -M 1499 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опция <span class="samp">-M</span> указывает программе максимальную длину последовательности, которую программа подаст на выход.</p>
<p><strong>Пример 3</strong>. Выбрать все последовательности, длина которых находится в промежутке [1500, 1600] п.н.</p>
<p class="samp-highl samp-vwide">seqkit seq -m 1500 -M 1600 ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<h4>Извлечение информации из заголовков последовательностей</h4>
<div id="select-head" class="pad-anchor"></div>
<p>Иногда требуется извлечь лишь данные fasta-заголовков, без самих последовательностей. Следующие примеры извлечения данных с помощью <span class="samp">seqkit seq</span> заголовков могут быть вам полезны:</p>
<p><strong>Пример 1</strong>. Выбрать все заголовки</p>
<p class="samp-highl samp-vwide">seqkit seq -n ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опция <span class="samp">-n</span> указывает программе подавать на выход только fasta-заголоки.</p>
<p><strong>Пример 2</strong>. Выбрать все идентификаторы последовательностей (части заголовков до первого пробела)</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz</p>
<p class="samp-comment">Опция <span class="samp">-i</span> указывает программе подавать на выход только идентификаторы последовательностей.</p>
<p><strong>Пример 3</strong>. Выбрать все коды доступа из базы данных RefSeq</p>
<p class="samp-highl samp-vwide">seqkit seq -ni ribogrove_{{ ribogrove_release_number }}_sequences.fasta.gz | sort | uniq</p>
<p class="samp-comment">Такая команда сработает успешно, лишь если в вашей системе установлены пррграммы <span class="samp">sort</span> и <span class="samp">uniq</span> (в системах Linux и Mac OS эти программы должны быть предустановлены).</p>
