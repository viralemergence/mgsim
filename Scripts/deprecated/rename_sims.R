rename_sims <- function(old_dir, new_dir, new_numbers) {
  successful <- old_dir %>% list.files() %>% gtools::mixedsort() %>%
    str_extract("[0-9]+") %>% as.numeric()
  new_paths <- new_numbers[successful] %>%
    map_chr(~file.path(new_dir, paste0("sample_", ., "_results.qs")))
  fs::file_move(old_dir %>% list.files(full.names = T) %>% gtools::mixedsort(),
                new_paths)
}
