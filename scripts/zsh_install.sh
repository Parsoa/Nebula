#!/bin/sh
# Build Zsh from sources on Ubuntu.
# From http://zsh.sourceforge.net/Arc/git.html and sources INSTALL file.

./Util/preconfig

# Options from Ubuntu Zsh package rules file (http://launchpad.net/ubuntu/+source/zsh)
./configure --prefix=$HOME/local \
            --enable-maildir-support \
            # --enable-function-subdirs \
            # --enable-site-fndir=$HOME/local/share/zsh/site-functions \
            # --enable-fndir=$HOME/share/zsh/functions \
            --with-tcsetpgrp \
            --with-term-lib="ncursesw" \
            --enable-cap \
            --enable-pcre \
            --enable-readnullcmd=pager \
            --enable-custom-patchlevel=Debian \
            LDFLAGS="-Wl,--as-needed -g"

make

make check

sudo make install

sudo make install.info

./configure --prefix=$HOME/local --enable-maildir-support --with-tcsetpgrp --with-term-lib="ncurses" --enable-cap --enable-pcre --enable-readnullcmd=pager --enable-custom-patchlevel=Debian LDFLAGS="-Wl,--as-needed -g" CFLAGS="-I/home/pkhorsand/local/include"