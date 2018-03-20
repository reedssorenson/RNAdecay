all: general_dExp_1sse.so general_Exp_1sse.so general_dExp_2sse.so general_Exp_2sse.so general_dExp_3sse.so general_Exp_3sse.so general_dExp_4sse.so general_Exp_4sse.so

general_dExp_1sse.so: general_dExp_1sse.cpp
	Rscript --vanilla -e "TMB::compile('general_dExp_1sse.cpp',safebounds=FALSE,safeunload=FALSE)"

general_Exp_1sse.so: general_Exp_1sse.cpp
	Rscript --vanilla -e "TMB::compile('general_Exp_1sse.cpp',safebounds=FALSE,safeunload=FALSE)"

general_dExp_2sse.so: general_dExp_2sse.cpp
	Rscript --vanilla -e "TMB::compile('general_dExp_2sse.cpp',safebounds=FALSE,safeunload=FALSE)"

general_Exp_2sse.so: general_Exp_2sse.cpp
	Rscript --vanilla -e "TMB::compile('general_Exp_2sse.cpp',safebounds=FALSE,safeunload=FALSE)"

general_dExp_3sse.so: general_dExp_3sse.cpp
	Rscript --vanilla -e "TMB::compile('general_dExp_3sse.cpp',safebounds=FALSE,safeunload=FALSE)"

general_Exp_3sse.so: general_Exp_3sse.cpp
	Rscript --vanilla -e "TMB::compile('general_Exp_3sse.cpp',safebounds=FALSE,safeunload=FALSE)"

general_dExp_4sse.so: general_dExp_4sse.cpp
	Rscript --vanilla -e "TMB::compile('general_dExp_4sse.cpp',safebounds=FALSE,safeunload=FALSE)"

general_Exp_4sse.so: general_Exp_4sse.cpp
	Rscript --vanilla -e "TMB::compile('general_Exp_4sse.cpp',safebounds=FALSE,safeunload=FALSE)"

clean:
	rm -rf *o
