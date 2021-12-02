reduction=c(0.3)


gr1c=0.15
gr2c=0.5
gr3c=0.23

gr1e=gr1c*(1-reduction)
gr2e=gr2c*(1-reduction)
gr3e=gr3c*(1-reduction)


ov32c=0.5*gr3c
ov32e=0.5*gr3e
ov123c=0.75*gr1c
ov123e=0.75*gr1e

gr23c=ov32c+gr2c-ov32c+gr3c-ov32c
gr23e=ov32e+gr2e-ov32e+gr3e-ov32e

gr123c=ov123c+gr1c-ov123c+gr23c-ov123c
gr123e=ov123e+gr1e-ov123e+gr23e-ov123e