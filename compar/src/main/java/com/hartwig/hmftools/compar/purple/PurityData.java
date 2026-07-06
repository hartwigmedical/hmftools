package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.common.CategoryType.PURITY;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.ComparableItem;

public class PurityData implements ComparableItem
{
    public final PurityContext Purity;

    protected static final String FLD_PURITY = "Purity";
    protected static final String FLD_PLOIDY = "Ploidy";
    protected static final String FLD_CONTAMINATION = "Contamination";
    protected static final String FLD_TMB = "TmbPerMb";
    protected static final String FLD_MS_INDELS = "MsIndelsPerMb";
    protected static final String FLD_TML = "Tml";
    protected static final String FLD_CN_SEGS = "CopyNumberSegments";
    protected static final String FLD_UNS_CN_SEGS = "UnsupportedCopyNumberSegments";
    protected static final String FLD_SV_TMB = "SvTmb";
    protected static final String FLD_QC_STATUS = "QcStatus";
    protected static final String FLD_GENDER = "Gender";
    protected static final String FLD_GERM_ABS = "GermlineAberrations";
    protected static final String FLD_FIT_METHOD = "FitMethod";
    protected static final String FLD_MS_STATUS = "MsStatus";
    protected static final String FLD_TMB_STATUS = "TmbStatus";
    protected static final String FLD_TML_STATUS = "TmlStatus";
    protected static final String FLD_TINC_LEVEL = "TincLevel";

    public PurityData(final PurityContext purityContext)
    {
        Purity = purityContext;
    }

    public CategoryType category() {
        return PURITY;
    }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }
}
