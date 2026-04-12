## Configuration

1. Rename  `config_template.json` to `config.json`.
2. Update `config.json` with your site-specific CLIF and output paths.
3. The current project expects, at minimum:
   - `site_name`
   - `clif_dir`
   - `output_dir`
   - `analysis_dir`
4. Optional exposome fields let you point the air-pollution workflow at county-level PM2.5, NO2, SVI, and ACS inputs.

Note: the `.gitignore` file in this directory ensures that the information in the config file is not pushed to github remote repository. 
