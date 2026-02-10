### TO DO

We can remove:

The function at ~202 (_update_info_dp) is called by format_for_singer, and that path is triggered when needs_singer_reformat(...) detects AD-style genotypes in the .clean records. That’s wired in flush_group (see split.py around the singer_reformat block).

If you no longer use the SINGER reformat path at all, we can safely remove:

needs_singer_reformat
format_for_singer
_update_info_dp
the singer_reformat logic in flush_group
Do you want me to remove that path entirely?



