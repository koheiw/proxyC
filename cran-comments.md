# Submission notes

## Purpose

Update documentation, add a new feature and fix a bug

## Test environments

* local: KDE neon User Edition 5.13, R 3.5.3
* local: Windows 10, R 3.5.3
* R-hub: Fedora Linux, R-devel, clang, gfortran (fedora-clang-devel)
* R-hub: Ubuntu Linux 16.04 LTS, R-devel, GCC (ubuntu-gcc-devel)
* R-hub: Ubuntu Linux 16.04 LTS, R-release, GCC (ubuntu-gcc-release)
* R-hub: Windows Server 2008 R2 SP1, R-devel, 32/64 bit (windows-x86_64-devel)
* R-hub: Windows Server 2008 R2 SP1, R-release, 32/64 bit (windows-x86_64-release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

Changes do not affect in existing users


```
> revdepcheck::revdep_reset()
> revdepcheck::revdep_check(timeout = as.difftime(600, units = "mins"), num_workers = 30)

todo: fill in

```
