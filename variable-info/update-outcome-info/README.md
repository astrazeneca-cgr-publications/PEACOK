
# Updating outcomes-info.tsv with new fields


UK Biobank describes its fields in a data dictionary on its website. When new fields are added we update the `../outcome-info.tsv` file used by PHESANT, using the updated UK Biobank data dictionary. 

We simplified the steps in original PHESANT package to work with the newer version of  UKB data dictionary file.

The steps we take to do this are:

1. Run pipeline/updateVariableInfor.r to download the latest UKB Data_Dictionary_Showcase.csv, make updated info file called `outcome-info-new.tsv`. 

This also makes a file `new-field-list.tsv` which lists the new fields added to `outcome-info-new.tsv` as well as any possible fields deleted. A copy of the downloaded Data_Dictionary_Showcase.csv will be kept for reference. 

2. Review fields and manually update PHESANT properties.

Fields with X in the PHESANT-specific columns need to be manually reviewed and values set in these columns as appropriate.

Fields with capitalized names in `../outcome-info.tsv` are additional PHESANT fields, used to process fields appropriately when running PHESANT.

And then move `outcome-info-new.tsv` to `../outcome-info.tsv`.



