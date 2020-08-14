package com.hartwig.hmftools.paddle.gene


data class GeneDriverRate(val gene: Gene, val nMissense: Double, val nNonsense: Double, val nSplice: Double, val nIndel: Double)