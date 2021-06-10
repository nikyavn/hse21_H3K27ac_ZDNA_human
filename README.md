# hse21_H3K27ac_ZDNA_human
# Проект по майнору "Биоинформатика", ВШЭ 2021
## hse21_H3K9me3_G4_human
*Усенкова Мария, 3 группа*
### Исходные данные

| Организм | Структура ДНК | Гистоновая метка | Тип клеток | Метка 1 | Метка 2 |
| -------- | ------------- | ---------------- | ---------- | ------- | ------- |
| Human (hg19) | G4_seq_Li_K | H3K9me3 | H1 | [ENCFF587TWB](https://www.encodeproject.org/files/ENCFF046DTX/) | [ENCFF697NMG](https://www.encodeproject.org/files/ENCFF697NMG/) |

Сохраненная сессия в UCSC GenomeBrowser: http://genome.ucsc.edu/s/mausenkova/hse21_H3K9me3_G4_human

### Анализ пиков гистоновой метки
Для работы были скачаны на кластер архивы с .bed-файлами с данными. При распаковке архивов были оставлены только первые 5 столбцов данных:
```bash
wget https://www.encodeproject.org/files/ENCFF697NMG/@@download/ENCFF697NMG.bed.gz
wget https://www.encodeproject.org/files/ENCFF587TWB/@@download/ENCFF587TWB.bed.gz
zcat ENCFF587TWB.bed.gz | cut -f1-5>H3K9me3_H1.ENCFF587TWB.hg19.bed
zcat ENCFF697NMG.bed.gz | cut -f1-5>H3K9me3_H1.ENCFF697NMG.hg19.bed
```
*Так как изначально скачанные .bed файлы были версии hg19, утилиту liftover я не использовала.*

Полученные файлы с помощью программы WinSCP были перенесены на ПК для дальнейшей работы.
##### Построение гистограмм длин участков
С помощью [скрипта](src/len_hist.R) на R были получены гистограммы длин участков для каждого эксперимента. 

Результаты:

![len_hist.ENCFF697NMG.hg19](images/png/len_hist.H3K9me3_H1.ENCFF697NMG.hg19.png)

![len_hist.ENCFF587TWB.hg19](images/png/len_hist.H3K9me3_H1.ENCFF587TWB.hg19.png)

##### Фильтрация пиков
С помощью [скрипта](/src/filtered.R) на R были отфильтрованы пики длиной более 5000. 

Результаты:

![filter_peaks.ENCFF697NMG.hg19.filtered.hist](images/png/filter_peaks.H3K9me3_H1.ENCFF697NMG.hg19.filtered.hist.png)

![filter_peaks.ENCFF587TWB.hg19.filtered.hist](images/png/filter_peaks.H3K9me3_H1.ENCFF587TWB.hg19.filtered.hist.png)

##### Расположение пиков

С помощью [скрипта](src/ChipSeeker.R) на R были построены графики расположения пиков гистоновых меток относительно аннотированных генов. 

Результаты:

###### chip_seeker.ENCFF697NMG.hg19.filtered.plotAnnoPie
![chip_seeker.ENCFF697NMG.hg19.filtered.plotAnnoPie](images/chip_seeker.H3K9me3_H1.ENCFF697NMG.hg19.filtered.plotAnnoPie.png)
###### chip_seeker.ENCFF587TWB.hg19.filtered.plotAnnoPie
![chip_seeker.ENCFF587TWB.hg19.filtered.plotAnnoPie](images/chip_seeker.H3K9me3_H1.ENCFF587TWB.hg19.filtered.plotAnnoPie.png)

##### Объединение файлов

Отсортированные файлы были загружены на кластер, отсортированы и объединены с помощью bedtools:

```bash
 cat  *.filtered.bed  |   sort -k1,1 -k2,2n   |   bedtools merge   >  H3K9me3_H1.merge.hg19.bed
```

Затем с помощью winSCP полученный файл был перенесен на ПК для дальнейшей работы.

##### Визуализация

С помощью [Genome Browser](http://genome.ucsc.edu/s/mausenkova/hse21_H3K9me3_human) были визуализированы полученные исходные наборы ChIP-seq пиков и их объединение:

```
track visibility=dense name="ENCFF587TWB"  description="H3K9me3_H1.ENCFF587TWB.hg19.filtered.bed"

track visibility=dense name="ENCFF697NMG"  description="H3K9me3_H1.ENCFF697NMG.hg19.filtered.bed"

track visibility=dense name="ChIP_merge"  color=50,50,200   description="H3K9me3_H1.merge.hg19.bed"

```

![GenomeBrowser1](images/png/GenomeBrowser1.png)

Объединение покрывает все наборы.

#### Анализ участков вторичной структуры ДНК

На кластер были скачены архивы с .bed-файлами с данными вторичной структуры ДНК:

```bash
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3003nnn/GSM3003539/suppl/GSM3003539_Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed.gz

wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3003nnn/GSM3003539/suppl/GSM3003539_Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed.gz
```

Распаковываем и удаляем не нужные для работы столбцы. Объединены в один файл с помощью bedtools:

```bash
zcat GSM3003539_Homo_all_w15_th-1_minus.hits.max.K.w50.25.bed.gz | cut -f1-5 > GSM3003539_minus.bed
zcat GSM3003539_Homo_all_w15_th-1_plus.hits.max.K.w50.25.bed.gz | cut -f1-5 > GSM3003539_plus.bed

cat GSM3003539_*.bed | sort -k1,1 -k2,2n | bedtools merge > GSM3003539.merged.bed 
```
Полученный файл с помощью программы WinSCP были перенесены на ПК для дальнейшей работы.

Далее с помощью [скрипта](src/len_hist.R) на R была получена гистограмма длин участков. 
Результаты:

![len_hist.GSM3003539.merged](images/png/len_hist.GSM3003539.merged.png)

Также с помощью [скрипта](src/ChipSeeker.R) на R был построен график расположения пиков относительно аннотированных генов.

Результаты:
###### chip_seeker.GSM3003539.merged.plotAnnoPie
![chip_seeker.GSM3003539.merged.plotAnnoPie](images/chip_seeker.GSM3003539.merged.plotAnnoPie.png)

#### Анализ пересечений гистоновой метки и структуры ДНК

С помощью bedtools были найдены пересечения гистоновой метки со структурами ДНК:
```bash
bedtools intersect  -a GSM3003539.merged.bed   -b  H3K9me3_H1.merge.hg19.bed  >  H3K9me3_H1.intersect_with_G4.bed
```

Полученный файл с помощью программы WinSCP были перенесены на ПК для дальнейшей работы.

Далее с помощью [скрипта](src/len_hist.R) на R была получена гистограмма длин участков. 
Результаты:

![H3K9me3_H1.intersect_with_G4](images/png/len_hist.H3K9me3_H1.intersect_with_G4.png)

Также с помощью [скрипта](src/ChipSeeker.R) на R был построен график расположения пиков относительно аннотированных генов.

Результаты:
###### chip_seeker.H3K9me3_H1.intersect_with_G4.plotAnnoPie
![chip_seeker.H3K9me3_H1.intersect_with_G4.plotAnnoPie](images/chip_seeker.H3K9me3_H1.intersect_with_G4.plotAnnoPie.png)

С помощью [Genome Browser](http://genome.ucsc.edu/s/mausenkova/hse21_H3K9me3_G4_human) были визуализированы полученные участки:

```
track visibility=dense name="ENCFF587TWB"  description="H3K9me3_H1.ENCFF587TWB.hg19.filtered.bed"

track visibility=dense name="ENCFF697NMG"  description="H3K9me3_H1.ENCFF697NMG.hg19.filtered.bed"

track visibility=dense name="ChIP_merge"  color=50,50,200   description="H3K9me3_H1.merge.hg19.bed"

track visibility=dense name="G4"  color=0,200,0  description="G4_Li_K"

track visibility=dense name="intersect_with_G4"  color=200,0,0  description="H3K9me3_H1.intersect_with_G4.bed"

```

Скриншот иллюстрирует пересечения между гистоновой меткой и структурой ДНК:

![GenomeBrowser2](images/png/GenomeBrowser2.png)

Например:

| Позиция | Координаты |
| ------- | ---------- |
| 1 | chr1:13790959-13791063 |
| 2 | chr1:13791788-13791843 |

Далее с помощью [скрипта](src/ChIPpeakAnno.R) на R полученные пересечения были ассоциированы с ближайшими генами. Было проассоциировано 490 пиков, 354 уникальных гена.

С помощью [Panther](http://pantherdb.org/) был проведён GO-анализ для полученных уникальных генов. В [файле](data/pantherdb_GO_analysis.txt) представлен результат анализа. 
Далее приведены значимые категории(c минимальными значениями FDR):
![pantherdb_result](images/png/pantherdb_result.png)

