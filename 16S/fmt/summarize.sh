#!/usr/bin/env bash

qiime demux summarize \
  --i-data fmt-tutorial-demux-1.qza \
  --o-visualization demux-summary-1.qzv
qiime demux summarize \
  --i-data fmt-tutorial-demux-2.qza \
  --o-visualization demux-summary-2.qzv
