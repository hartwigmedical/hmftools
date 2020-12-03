package com.hartwig.hmftools.gripss

enum class SvType(val includeInVcf: Boolean) {
    DEL(true), INS(true), DUP(true), INV(true), SGL(false), BND(true);
}