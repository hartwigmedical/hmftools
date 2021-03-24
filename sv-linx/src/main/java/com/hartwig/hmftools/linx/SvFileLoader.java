package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.valueNotNull;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFile;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.linx.types.GermlineFilter;
import com.hartwig.hmftools.patientdb.dao.DatabaseUtil;

public class SvFileLoader
{
    public static final String VCF_FILE = "sv_vcf";

    public static List<StructuralVariantData> loadSvDataFromVcf(final String vcfFile)
    {
        final List<StructuralVariantData> svDataList = Lists.newArrayList();

        try
        {
            final List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter());
            final List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

            // generate a unique ID for each SV record
            int svId = 0;

            for (EnrichedStructuralVariant var : enrichedVariants)
            {
                svDataList.add(convertSvData(var, svId++));
            }
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to load SVs from VCF: {}", e.toString());
        }

        LNX_LOGGER.info("loaded {} SV data records from VCF file: {}", svDataList.size(), vcfFile);

        return svDataList;
    }

    public static List<StructuralVariantData> loadSvDataFromGermlineVcf(final String vcfFile)
    {
        final List<StructuralVariantData> svDataList = Lists.newArrayList();

        try
        {
            final List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfFile, new GermlineFilter(true));

            int svId = 0;

            for (StructuralVariant var : variants)
            {
                svDataList.add(convertGermlineSvData(var, svId++));
            }
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to load SVs from VCF: {}", e.toString());
        }

        LNX_LOGGER.info("loaded {} germline SV data records from VCF file: {}", svDataList.size(), vcfFile);

        return svDataList;
    }

    public static List<StructuralVariantData> loadSvDataFromSvFile(final String sampleId, final String svDataPath)
    {
        try
        {
            final String svDataFile = StructuralVariantFile.generateFilename(svDataPath, sampleId);
            return StructuralVariantFile.read(svDataFile);
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to load SV data: {}", e.toString());
            return Lists.newArrayList();
        }
    }

    public static StructuralVariantData convertGermlineSvData(final StructuralVariant var, int svId)
    {
        return ImmutableStructuralVariantData.builder()
                .id(svId)
                .startChromosome(var.chromosome(true))
                .endChromosome(var.end() == null ? "0" : var.chromosome(false))
                .startPosition(var.position(true).intValue())
                .endPosition(var.end() == null ? -1 : var.position(false).intValue())
                .startOrientation(var.orientation(true))
                .endOrientation(var.end() == null ? (byte) 0 : var.orientation(false))
                .startHomologySequence(var.start().homology())
                .endHomologySequence(var.end() == null ? "" : var.end().homology())
                .junctionCopyNumber(1)
                .startAF(DatabaseUtil.valueNotNull(var.start().alleleFrequency()))
                .endAF(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().alleleFrequency()))
                .adjustedStartAF(DatabaseUtil.valueNotNull(var.start().alleleFrequency()))
                .adjustedEndAF(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().alleleFrequency()))
                .adjustedStartCopyNumber(DatabaseUtil.valueNotNull(1))
                .adjustedEndCopyNumber(var.end() == null ? 0 : 1)
                .adjustedStartCopyNumberChange(1)
                .adjustedEndCopyNumberChange(var.end() == null ? 0 : 1)
                .insertSequence(var.insertSequence())
                .type(var.type())
                .filter(var.filter())
                .imprecise(var.imprecise())
                .qualityScore(DatabaseUtil.valueNotNull(var.qualityScore()))
                .event(valueNotNull(var.event()))
                .startTumorVariantFragmentCount(DatabaseUtil.valueNotNull(var.start().tumorVariantFragmentCount()))
                .startTumorReferenceFragmentCount(DatabaseUtil.valueNotNull(var.start().tumorReferenceFragmentCount()))
                .startNormalVariantFragmentCount(DatabaseUtil.valueNotNull(var.start().normalVariantFragmentCount()))
                .startNormalReferenceFragmentCount(DatabaseUtil.valueNotNull(var.start().normalReferenceFragmentCount()))
                .endTumorVariantFragmentCount(0)
                .endTumorReferenceFragmentCount(0)
                .endNormalVariantFragmentCount(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().normalVariantFragmentCount()))
                .endNormalReferenceFragmentCount(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().normalReferenceFragmentCount()))
                .startIntervalOffsetStart(DatabaseUtil.valueNotNull(var.start().startOffset()))
                .startIntervalOffsetEnd(DatabaseUtil.valueNotNull(var.start().endOffset()))
                .endIntervalOffsetStart(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().startOffset()))
                .endIntervalOffsetEnd(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().endOffset()))
                .inexactHomologyOffsetStart(DatabaseUtil.valueNotNull(var.start().inexactHomologyOffsetStart()))
                .inexactHomologyOffsetEnd(DatabaseUtil.valueNotNull(var.start().inexactHomologyOffsetEnd()))
                .startLinkedBy(valueNotNull(var.startLinkedBy()))
                .endLinkedBy(valueNotNull(var.endLinkedBy()))
                .vcfId(valueNotNull(var.id()))
                .startRefContext("") // getValueNotNull(var.start().refGenomeContext())
                .endRefContext(var.end() == null ? "" : "") // getValueNotNull(var.end().refGenomeContext())
                .recovered(var.recovered())
                .recoveryMethod((valueNotNull(var.recoveryMethod())))
                .recoveryFilter(valueNotNull(var.recoveryFilter()))
                .insertSequenceAlignments(valueNotNull(var.insertSequenceAlignments()))
                .insertSequenceRepeatClass(valueNotNull(var.insertSequenceRepeatClass()))
                .insertSequenceRepeatType(valueNotNull(var.insertSequenceRepeatType()))
                .insertSequenceRepeatOrientation(DatabaseUtil.valueNotNull(var.insertSequenceRepeatOrientation()))
                .insertSequenceRepeatCoverage(DatabaseUtil.valueNotNull(var.insertSequenceRepeatCoverage()))
                .startAnchoringSupportDistance(var.start().anchoringSupportDistance())
                .endAnchoringSupportDistance(var.end() == null ? 0 : var.end().anchoringSupportDistance())
                .build();
    }
}
