# Comparing enhancer activity between KO and WT cell lines
Steps to detect changes in enhancer activity in the KO line compared to the WT line are as follows:

1. Filter the original master file and keep only those regions which had greater than **50 filtered reads** assigned to them per replicate of the Input STARR-seq library.
    
    src code: root/src/0_create_filtered_master_list.py

2. Break the filtered master list into overlapping windows of size 500 and stride 50.

    src code: root/src/1_create_overlapping_windows.py

3. Calculate filtered reads assigned to the overlapping windows per replicate of all libraries.

    src code: root/src/2_get_window_depth.py

4. Compare the filtered reads assigned to the overlapping windows between KO and WT libraries using DESeq2

    src code: root/src/
