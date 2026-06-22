package com.hartwig.hmftools.orange.report.pdfdata;

import java.util.List;

import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.purple.PurpleChrArmCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;

import org.jetbrains.annotations.Nullable;

public class SomaticFindingsData
{
    public final boolean hasPurpleFail;
    public final boolean tumorOnlyMode;
    public final boolean hasRna;
    public final boolean hasRnaSample;
    public final List<PurpleVariant> somaticVariants;
    public final List<PurpleGainDeletion> somaticGainsDels;
    public final List<PurpleChrArmCopyNumber> armCopyNumberAbberations;
    public final List<LinxFusion> fusions;
    public final List<LinxBreakend> somaticBreakends;
    @Nullable
    public final VirusInterpreterData virusInterpreter;
    @Nullable
    public final List<SignatureAllocation> sigAllocations;
    public final List<String> linxDriverPlotPaths;

    public SomaticFindingsData(
            final boolean hasPurpleFail,
            final boolean tumorOnlyMode,
            final boolean hasRna,
            final boolean hasRnaSample,
            final List<PurpleVariant> somaticVariants,
            final List<PurpleGainDeletion> somaticGainsDels,
            final List<PurpleChrArmCopyNumber> armCopyNumberAbberations,
            final List<LinxFusion> fusions,
            final List<LinxBreakend> somaticBreakends,
            @Nullable final VirusInterpreterData virusInterpreter,
            @Nullable final List<SignatureAllocation> sigAllocations,
            final List<String> linxDriverPlotPaths)
    {
        this.hasPurpleFail = hasPurpleFail;
        this.tumorOnlyMode = tumorOnlyMode;
        this.hasRna = hasRna;
        this.hasRnaSample = hasRnaSample;
        this.somaticVariants = somaticVariants;
        this.somaticGainsDels = somaticGainsDels;
        this.armCopyNumberAbberations = armCopyNumberAbberations;
        this.fusions = fusions;
        this.somaticBreakends = somaticBreakends;
        this.virusInterpreter = virusInterpreter;
        this.sigAllocations = sigAllocations;
        this.linxDriverPlotPaths = linxDriverPlotPaths;
    }
}
