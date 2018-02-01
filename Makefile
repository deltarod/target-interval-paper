FPOP.gridsearch.RData: FPOP.gridsearch.R
	R --no-save < $<
figure-FPOP-PDPA-compare.png: figure-FPOP-PDPA-compare.R FPOP.models.RData PDPA.models.RData
	R --no-save < $<
FPOP.models.RData: FPOP.models.R
	R --no-save < $<
PDPA.models.RData: PDPA.models.R
	R --no-save < $<
