Check of the repeat track by the main session, 2026-11-19
==========================================================

Why: the author asked for an explanation of the note on the track that the draft calls
"transposable elements" (03_input_files/genome_features/cgigas_uk_roslin_v1_rm.te.bed.gz)
and whether the list should be curated. The numbers of that note came from one agent of the
code pass (group B) and had not been reviewed. This is an independent recount. It changes
no pipeline file and no output.

Inputs
- NCBI's RepeatMasker output for the assembly:
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/806/645/GCF_902806645.1_cgigas_uk_roslin_v1/GCF_902806645.1_cgigas_uk_roslin_v1_rm.out.gz
  (md5 of the downloaded file e3c515b1edfb9941d9b60012c395aa75; not held in this repository)
- united_cpg_legacy.csv.gz, dml_ploidy_legacy.csv, dml_heat_legacy.csv, dml_shared_legacy.csv of
  04_output/05_differential_methylation/tables, as committed (taken with git show HEAD:...).

What was checked, and the result
1. Columns 5 to 7 of the RepeatMasker output, without its three header lines, are the committed
   track line for line (diff of the two: no difference; 781,312 records).
2. Records and bases by class, counted with awk (class = column 11 before the slash):
   DNA 380,683 records, 116,812,111 bp; Simple_repeat 165,597, 8,116,222; RC (all RC/Helitron)
   151,061, 84,257,716; LINE 29,494, 13,491,361; Low_complexity 21,815, 1,050,398; LTR 19,774,
   13,555,792; Unknown 9,973, 3,641,919; Satellite 1,461, 363,422; SINE 1,391, 242,275; rRNA 63,
   12,887. Total 781,312 records, 241,544,103 bp (sum over records). Simple repeats, low
   complexity, satellites and rRNA: 188,936 records (24.18%), 3.95% of the bases.
   190,162 of the DNA records carry no family name (class/family is "DNA" alone).
3. The coordinates are 1-based and closed. The first record, (C)n at 4109 to 4128 of NC_047559.1,
   matches positions 1 to 20 of the repeat; in the reference, bases 4109 to 4128 (1-based) are
   twenty C, preceded by "a" and followed by "a". Read as a BED interval the record covers
   nineteen of them. BED: chromStart is 0-based and chromEnd is not included
   (https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
4. Overlap of the analysed CpGs (legacy mode) with the track, by foverlaps of data.table, a CpG
   being the two bases start + 1 and start + 2 (te_classes_check.R; output below). Agrees with
   the numbers of group B (te_readings.R, verify_tracks.R of the code pass):
   read as BED, all classes 109,966 CpGs, 194 ploidy DML, 463 heat DML, 28 loci in both lists;
   coordinates as written 110,067, 194, 465, 28;
   classes DNA, RC, LTR, LINE and SINE only: 106,269 (as BED) or 106,342 (as written) CpGs,
   189 ploidy DML, 455 or 456 heat DML, 28 in both lists.
   New here: of the 3,697 CpGs outside those five classes (as BED), 2,158 lie in a record of
   class Unknown and 1,539 in simple repeats, low complexity, satellites or rRNA alone; of the
   5 ploidy DML outside the five classes 3 are in Unknown records and 2 in the others; of the
   8 heat DML, 5 and 3.

How NCBI made the file (GCF_902806645.1_cgigas_uk_roslin_v1_rm.run, same folder at NCBI, read
2026-11-19):
----
The repeat files provided for this assembly were generated using RepeatMasker.
  Smit, AFA, Hubley, R & Green, P.,
  RepeatMasker Open-3.0.
  1996-2010 <http://www.repeatmasker.org>.

VERSION:
RepeatMasker version open-4.0.8 , sensitive mode
run with blastp version 2.0MP-WashU [01-Jan-2006] [linux24-i786-ILP32F64 2006-01-02T05:13:21]
RepeatMasker Combined Database: Dfam_Consensus-20181026, RepBase-20181026

PARAMETERS:
RepeatMasker -engine wublast -species 'crassostrea gigas' -s -no_is -cutoff 255 -frag 20000

REPEATS:
RepeatMasker Database: RepeatMaskerLib.embl
Version: RepeatMasker Combined Database: Dfam_Consensus-20181026, RepBase-20181026
Species: crassostrea gigas ( crassostrea gigas )
180 ancestral and ubiquitous sequence(s) with a total length of 63759 bp
798 crassostrea gigas specific repeats with a total length of 2041595 bp
0 lineage specific sequence(s) with a total length of 0 bp
--------------------------------------------------------------------------------
----
RepeatMasker "screens DNA sequences for interspersed repeats and low complexity DNA sequences"
(https://www.repeatmasker.org/webrepeatmaskerhelp.html); the run above has no option that
leaves simple repeats or low complexity out, so the output holds them.

Output of te_classes_check.R
----
===== reading: as BED (legacy)
            set n_loci all_classes   five five_or_unknown unknown_not_five
1: all analysed 480579      109966 106269          108427             2158
2:   ploidy DML   1038         194    189             192                3
3:     heat DML   2351         463    455             460                5
4:       shared    200          28     28              28                0
   only_simple_low_sat_rrna
1:                     1539
2:                        2
3:                        3
4:                        0
DML outside the five classes, by class combination:
         classes    set n
1: Simple_repeat ploidy 2
2:       Unknown ploidy 3
3: Simple_repeat   heat 2
4:       Unknown   heat 5
5:     Satellite   heat 1

===== reading: coordinates as written
            set n_loci all_classes   five five_or_unknown unknown_not_five
1: all analysed 480579      110067 106342          108502             2160
2:   ploidy DML   1038         194    189             192                3
3:     heat DML   2351         465    456             462                6
4:       shared    200          28     28              28                0
   only_simple_low_sat_rrna
1:                     1565
2:                        2
3:                        3
4:                        0
DML outside the five classes, by class combination:
         classes    set n
1: Simple_repeat ploidy 2
2:       Unknown ploidy 3
3: Simple_repeat   heat 2
4:       Unknown   heat 6
5:     Satellite   heat 1

mean record length by class (bp):
             class records mean_bp median_bp
 1:            DNA  380683   306.8       201
 2:  Simple_repeat  165597    49.0        40
 3:             RC  151061   557.8       272
 4:           LINE   29494   457.4       208
 5: Low_complexity   21815    48.2        42
 6:            LTR   19774   685.5       296
 7:        Unknown    9973   365.2       192
 8:      Satellite    1461   248.7       210
 9:           SINE    1391   174.2       194
10:           rRNA      63   204.6        63

DNA records without a family: 190162  Unknown: 9973
----
