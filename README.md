# EasyTargetscan

EasyTargetScan_0.1.py is a simple script that allows to use fasta files directly in with TargetScan to predict miRNA binding sites on a given sequence.

it requires Pearl and TargetScan
http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi

http://strawberryperl.com/

```
usage: EasyTargetScan_0.1.py [-h] -utr UTR -mir MIR [-output OUTPUT]
                             [-Taxon TAXON]

Allows to use TargetScan from FASTA files

optional arguments:
  -h, --help      show this help message and exit
  -utr UTR        a FASTA file conataining sequences to scan
  -mir MIR        a FASTA file conataining mature miRNA sequences
  -output OUTPUT  prefix of the output file
  -Taxon TAXON    mouse by default:10090
```

A simple scan using seed sequences can be done with the second python script scan_seq_with_seeds.py (without Pearl and TargetScan)
