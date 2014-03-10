Tests for InterFits / LedaFits
------------------------------

0) Run corr2uvfits.py - This converts test_lalc.LA and .LC into a uvfits file, which
   is used a "gold standard".
   
1) test_fitsidi.py - This tests that conversion from uvfits -> fitsidi works as intended.
   It generates a fitsidi file from the uvfits file, then generates a new fitsidi file to
   confirm that uvfits -> interfits -> fitsidi -> interfits -> fitsidi does not lose data.
   After this, the test_lalc.fitsidi file is considered the "gold standard"

2) test_fitsidi_lalc.py - This tests that reading the LA file directly into interfits creates
   valid flux data. 

3) test_lalc_direct.py - This tests that conversion from LA->interfits->fitsidi->interfits gives
   the same data as reading LA->interfits
   
4) test_dada_direct.py - This tests the reading of a dada file and conversion into fitsidi, as 
   done in test_lalc_direct