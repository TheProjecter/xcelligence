##----------------------------------------------------------------------------##
##
##	Makefile, created on Jul 13, 2009
##	Author: Jitao David Zhang <j.zhang@dkfz.de>
##	Description: makefile for building distributions etc.
##
##	The Makefile provides the following targets:
##                  
##      - make install  	calls R-devel CMD INSTALL
##	- make stable-install	calls R CMD INSTALL
##      - make check    	calls R-devel CMD check (with RUnit)
##	- make stable-check	calls R CMD check (with RUnit)
##      - make dist     	calls R-devel CMD build
##      - make stable-dist     	calls R CMD build
##			
##
##----------------------------------------------------------------------------##

R=R-devel
Rstable=R
PKG          := xcelligence
PKG_VERSION  := dist

PKG_ROOT_DIR := $(shell pwd)
PKG_DIST_ROOT_DIR := ../$(PKG).tmp
PKG_HIDDEN_FILES  := Makefile 

install: 
	@echo '====== Installing Package ======'
	@(cd ..; ${R} CMD INSTALL $(PKG))
	@echo '====== Installing finished ======'
	@echo ' '

stable-install: 
	@echo '====== Installing Package ======'
	@(cd ..; ${Rstable} CMD INSTALL $(PKG))
	@echo '====== Installing finished ======'
	@echo ' '

check:	
	@echo '====== Checking Package ======'
	@(export R_DEVELOP_MODE=TRUE; cd ..; ${R} CMD check $(PKG))
	@echo '====== Checking finished ======'
	@echo ' '

stable-check:	
	@echo '====== Checking Package ======'
	@(export R_DEVELOP_MODE=TRUE; cd ..; ${Rstable} CMD check $(PKG))
	@echo '====== Checking finished ======'
	@echo '

dist:	
	@echo '====== Building Distribution ======'
	cp -rp $(PKG_ROOT_DIR) $(PKG_DIST_ROOT_DIR)
	@(cd ..; $(RM) -r $(PKG_HIDDEN_FILES); ${R} CMD build $(PKG))
	$(RM) -r $(PKG_DIST_ROOT_DIR)
	@echo '====== Building finished ======'
	@echo ' '
	@echo '====== Checking Package ======'
	@(export R_DEVELOP_MODE=TRUE; cd ..; ${R} CMD check $(PKG)_$(PKG_VERSION).tar.gz)
	@echo '====== Checking finished  ======'
	@echo ' '

stable-dist:	
	@echo '====== Building Distribution ======'
	cp -rp $(PKG_ROOT_DIR) $(PKG_DIST_ROOT_DIR)
	@(cd ..; $(RM) -r $(PKG_HIDDEN_FILES); ${Rstable} CMD build $(PKG))
	$(RM) -r $(PKG_DIST_ROOT_DIR)
	@echo '====== Building finished ======'
	@echo ' '
	@echo '====== Checking Package ======'
	@(export R_DEVELOP_MODE=TRUE; cd ..; ${Rstable} CMD check $(PKG)_$(PKG_VERSION).tar.gz)
	@echo '====== Checking finished  ======'
	@echo ' '