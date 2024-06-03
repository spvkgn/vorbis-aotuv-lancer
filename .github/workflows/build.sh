#!/bin/bash
set -eo pipefail

BUILD_DIR=$GITHUB_WORKSPACE/build

export CC="ccache gcc"
export CFLAGS="-O3 -fno-strict-aliasing"
export PKG_CONFIG_PATH="$BUILD_DIR/lib/pkgconfig"
export CCACHE_DIR="$GITHUB_WORKSPACE/.ccache"

case $ARCH in
  x86_64)  PLATFORM=x64 ;;
  x86)     PLATFORM=x86 ;;
  aarch64) PLATFORM=arm64 ;;
  armhf)   PLATFORM=arm ;;
  *)       PLATFORM=$ARCH ;;
esac

SUFFIX=-aotuv-lancer

mkdir -p "$BUILD_DIR"
cd $GITHUB_WORKSPACE

get_sources() {
  SOURCES_URL=http://downloads.xiph.org/releases/
  NAME=$1
  case $NAME in
    *vorbis*) SOURCES_URL+=vorbis ;;
    *ogg*)  SOURCES_URL+=ogg ;;
    *flac*) SOURCES_URL+=flac ;;
  esac
  echo "Download $NAME sources"
  wget -qO- $SOURCES_URL |\
    grep -Po "href=\"\K$NAME-(\d+\.)+\d+.*\.tar\.(gz|xz)(?=\")" | sort -V | tail -1 |\
    xargs -I{} wget -qO- $SOURCES_URL/{} | bsdtar -x
}

# build libogg
get_sources libogg
( cd libogg-*/
  autoreconf -fi && \
  ./configure --prefix=$BUILD_DIR \
    --disable-shared --enable-static \
    --disable-dependency-tracking && \
  make -j$(nproc) install )

# build FLAC
get_sources flac
( cd flac-*/
  autoreconf -fi && \
  ./configure --prefix=$BUILD_DIR \
    --disable-shared --enable-static \
    --disable-dependency-tracking \
    --disable-debug \
    --disable-oggtest \
    --disable-programs \
    --disable-examples \
    --disable-cpplibs \
    --disable-doxygen-docs \
    --with-ogg=$BUILD_DIR && \
  make -j$(nproc) install )

# build Vorbis
( cd vorbis-aotuv-lancer
  ./autogen.sh && \
  ./configure --prefix=$BUILD_DIR \
    --disable-shared --enable-static \
    --disable-dependency-tracking \
    --disable-maintainer-mode \
    --with-ogg=$BUILD_DIR && \
  make -j$(nproc) install )

# build vorbis-tools
get_sources vorbis-tools
( cd vorbis-tools-*/
  wget -qO- https://git.alpinelinux.org/aports/plain/community/vorbis-tools/CVE-2023-43361.patch | patch -p1 -i -
  autoreconf -fi && \
  if [[ $PLATFORM == arm ]]; then
    sed 's/cross_compiling=no/cross_compiling=yes/g' -i configure
  fi
  sed 's/@LDFLAGS@/@LDFLAGS@ -all-static/' -i Makefile.in && \
  CPPFLAGS="-I$BUILD_DIR/include" \
  LDFLAGS="-Wl,-static -static-libgcc -no-pie -L$BUILD_DIR/lib" \
  ./configure \
    --prefix=$BUILD_DIR \
    --program-suffix=$SUFFIX \
    --disable-dependency-tracking \
    --disable-maintainer-mode \
    --without-speex \
    --without-curl \
    --disable-ogg123 \
    --disable-oggdec \
    --disable-ogginfo \
    --disable-vcut \
    --disable-vorbiscomment \
    --with-ogg=$BUILD_DIR \
    --with-vorbis=$BUILD_DIR && \
  make -j $(nproc) install || exit 1 )

( cd $BUILD_DIR/bin ; strip oggenc* ; tar -cJvf $GITHUB_WORKSPACE/oggenc$SUFFIX-$PLATFORM.tar.xz oggenc* )

ccache --show-stats
