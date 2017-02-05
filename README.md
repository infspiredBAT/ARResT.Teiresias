# ARResT.TEIRESIAS
discovering new subsets of stereotyped antigen receptor sequences

</br>
citation: http://www.ncbi.nlm.nih.gov/pubmed/19759557</br>
code: we provide Perl scripts, IBM's TEIRESIAS (old version), NR-grep by G. Navarro, and demo data</br>
usage:
<pre>
teiresias.pl

--in		-> table: tab-delimited, NO header line, Unix ASCII - see demo.table
                 1st column:  required!  label, make sure they're unique
                 2nd column:  required!  amino acid junction sequence, with anchors
                 3rd column:  optional   IMGT-formatted V gene (and allele)
                 4th column:  optional   annotation

    other (optional) arguments - defaults are current standards for CLL (PMID:22415752)
--ide       -> minimum identity [0-1], default: 0.5
--sim       -> minimum similarity [0-1], default: 0.7
--lendiff   -> max length difference [0,1,2,3], default: 0
--offdiff   -> max offset difference [0,1,2,3], default: 0
--nophylo   -> ignore phylogenetic relationships as defined in 'phyloequiv'

--demo      -> run demo (~30sec in a multicore system)
--debug     -> do not erase intermediate files


    note / hints / tips
When possible, sequences in results table are uppercased with the pattern that clustered them.
You can find the same sequence in different clusters (i.e. ~fuzzy clustering)
 - these are used as linkers for higher-level clustering (L1, L2, etc),
   where linkage can actually be very loose so be careful with interpretation.
</pre>
