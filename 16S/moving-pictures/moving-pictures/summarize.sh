#!/usr/bin/env bash

qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv
