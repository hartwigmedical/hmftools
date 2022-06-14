package com.hartwig.hmftools.compar.purple;

import static com.hartwig.hmftools.compar.Category.PURITY;
import static com.hartwig.hmftools.compar.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.MatchLevel.REPORTABLE;
import static com.hartwig.hmftools.compar.MismatchType.VALUE;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.compar.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.DiffThresholds;
import com.hartwig.hmftools.compar.MatchLevel;
import com.hartwig.hmftools.compar.Mismatch;

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

    public PurityData(final PurityContext purityContext)
    {
        Purity = purityContext;
    }

    public Category category() {
        return PURITY;
    }

    @Override
    public String key()
    {
        return "";
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(String.format("%.2f", Purity.bestFit().purity()));
        values.add(String.format("%.2f", Purity.bestFit().ploidy()));
        values.add(String.format("%.4f", Purity.qc().contamination()));
        values.add(String.format("%.2f", Purity.tumorMutationalBurdenPerMb()));
        values.add(String.format("%d", Purity.tumorMutationalLoad()));
        values.add(String.format("%.4f", Purity.microsatelliteIndelsPerMb()));
        values.add(String.format("%d", Purity.svTumorMutationalBurden()));
        values.add(String.format("%d", Purity.qc().copyNumberSegments()));
        values.add(String.format("%d", Purity.qc().unsupportedCopyNumberSegments()));
        values.add(String.format("%s", qcStatus(Purity.qc().status())));
        values.add(String.format("%s", Purity.qc().cobaltGender()));
        values.add(String.format("%s", germlineAberrations(Purity.qc().germlineAberrations())));
        values.add(String.format("%s", Purity.method()));
        values.add(String.format("%s", Purity.microsatelliteStatus()));
        values.add(String.format("%s", Purity.tumorMutationalBurdenStatus()));
        values.add(String.format("%s", Purity.tumorMutationalLoadStatus()));

        return values;
    }

    @Override
    public boolean reportable() {
        return true;
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        // a single record for each sample
        return true;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final PurityData otherPurity = (PurityData) other;

        final List<String> diffs = Lists.newArrayList();

        checkDiff(diffs, FLD_PURITY, Purity.bestFit().purity(), otherPurity.Purity.bestFit().purity(), thresholds);
        checkDiff(diffs, FLD_PLOIDY, 1.0, otherPurity.Purity.bestFit().ploidy(), thresholds);

        checkDiff(diffs, FLD_CONTAMINATION, Purity.qc().contamination(), otherPurity.Purity.qc().contamination(), thresholds);

        checkDiff(diffs, FLD_TMB, Purity.tumorMutationalBurdenPerMb(), otherPurity.Purity.tumorMutationalBurdenPerMb(), thresholds);
        checkDiff(diffs, FLD_MS_INDELS, Purity.microsatelliteIndelsPerMb(), otherPurity.Purity.microsatelliteIndelsPerMb(), thresholds);
        checkDiff(diffs, FLD_TML, Purity.tumorMutationalLoad(), otherPurity.Purity.tumorMutationalLoad(), thresholds);

        checkDiff(diffs, FLD_CN_SEGS, Purity.qc().copyNumberSegments(), otherPurity.Purity.qc().copyNumberSegments(), thresholds);

        checkDiff(
                diffs, FLD_UNS_CN_SEGS,
                Purity.qc().unsupportedCopyNumberSegments(), otherPurity.Purity.qc().unsupportedCopyNumberSegments(),thresholds);

        checkDiff(diffs, FLD_SV_TMB, Purity.bestFit().purity(), otherPurity.Purity.bestFit().purity(), thresholds);

        checkDiff(diffs, FLD_QC_STATUS, qcStatus(Purity.qc().status()), qcStatus(otherPurity.Purity.qc().status()));

        checkDiff(diffs, FLD_GENDER, Purity.gender().toString(), otherPurity.Purity.gender().toString());

        checkDiff(
                diffs, FLD_GERM_ABS,
                germlineAberrations(Purity.qc().germlineAberrations()), germlineAberrations(otherPurity.Purity.qc().germlineAberrations()));

        checkDiff(diffs, FLD_FIT_METHOD, Purity.method().toString(), otherPurity.Purity.method().toString());

        checkDiff(diffs, FLD_MS_STATUS, Purity.microsatelliteStatus().toString(), otherPurity.Purity.microsatelliteStatus().toString());

        checkDiff(
                diffs, FLD_TMB_STATUS,
                Purity.tumorMutationalBurdenStatus().toString(), otherPurity.Purity.tumorMutationalBurdenStatus().toString());

        checkDiff(
                diffs, FLD_TML_STATUS,
                Purity.tumorMutationalLoadStatus().toString(), otherPurity.Purity.tumorMutationalLoadStatus().toString());

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }

    private static String germlineAberrations(final Set<GermlineAberration> aberrations)
    {
        StringJoiner sj = new StringJoiner(";");
        aberrations.forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public static String qcStatus(final Set<PurpleQCStatus> status)
    {
        StringJoiner sj = new StringJoiner(";");
        status.forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }
}
