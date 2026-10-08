0.5
* Add POLCA and NextPolish as polishing steps. Generalize polishing beyond pilon [Done]
* add FCS-GX support (w/ or w/o memory mapped approach) [Done]
* support masurca and nextdenovo assemblers [Done]

* consider flye or canu support

## Known issues

* **taxpurge memory**: `fastmultigather` loads the sketch collections into memory (GTDB rs226
  representatives + GenBank fungi together: ~1.1 billion hashes). The branchwater RocksDB indexes
  (`sourmash scripts index`, or the published `gtdb-rs226-reps.k31.rocksdb`) are searched on disk
  and would cut memory a lot; consider supporting them in `taxpurge --sourdb`.
