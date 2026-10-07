0.5
* Add POLCA and NextPolish as polishing steps. Generalize polishing beyond pilon [Done]
* add FCS-GX support (w/ or w/o memory mapped approach) [Done]
* support masurca and nextdenovo assemblers [Done]

* consider flye or canu support

## Known issues

* **taxpurge GTDB databases** (`aaftf/resources.py`, `sm_gtdb` / `sm_gtdbrep`): the URLs download
  sourmash `.dna.zip` signature collections, but the files are saved as `.lca.json.gz` and
  `taxpurge` runs `sourmash lca classify --db`, which needs an LCA database. So
  `--sourdb_type gtdb` / `gtdbrep` likely fail; the default `gbk` database is unaffected. Fix by
  pointing these entries at LCA databases, or by switching taxpurge to a classifier that accepts
  zip collections (e.g. `sourmash gather` + `sourmash tax`).
