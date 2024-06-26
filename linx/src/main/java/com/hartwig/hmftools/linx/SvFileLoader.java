package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.sv.SvVcfTags.INFERRED;
import static com.hartwig.hmftools.common.sv.SvVcfTags.PON_COUNT;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_CN_CHANGE;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_JUNCTION_COPY_NUMBER;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.gridss.GridssSvFactory;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.linx.germline.GermlineFilter;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.logging.log4j.util.Strings;

import htsjdk.variant.variantcontext.VariantContext;

public final class SvFileLoader
{
    public static List<StructuralVariantData> loadVariantsFromVcf(final LinxConfig config, final String sampleId)
    {
        String vcfFile = convertWildcardSamplePath(config.SvVcfFile, sampleId);

        if(config.IsGermline)
            return loadGermlineVariantsFromVcf(vcfFile, sampleId);
        else
            return loadSomaticVariantsFromVcf(vcfFile);
    }

    private static List<StructuralVariantData> loadSomaticVariantsFromVcf(final String vcfFile)
    {
        List<StructuralVariantData> svDataList = Lists.newArrayList();

        try
        {
            List<StructuralVariant> variants = StructuralVariantFileLoader.fromGridssFile(vcfFile, new AlwaysPassFilter());
            // List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter());
            List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

            // generate a unique ID for each SV record
            int svId = 0;

            for(EnrichedStructuralVariant var : enrichedVariants)
            {
                svDataList.add(convertSvData(var, svId++));
            }

            LNX_LOGGER.info("loaded {} SV data records from VCF file: {}", svDataList.size(), vcfFile);
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to load SVs from VCF: {}", e.toString());
        }

        return svDataList;
    }

    private static List<StructuralVariantData> loadGermlineVariantsFromVcf(final String vcfFile, final String sampleId)
    {
        // StructuralVariantFactory svFactory = StructuralVariantFactory.build(new GermlineFilter());
        GridssSvFactory svFactory = GridssSvFactory.build(new GermlineFilter());

        VcfFileReader vcfReader = new VcfFileReader(vcfFile);

        if(!vcfReader.fileValid())
            return Collections.emptyList();

        // NOTE: if linx (and Purple) are run in germline-only mode, then consider setting the reference ordinal only here,
        // but depends on the downstream interpretation of the Linx output
        List<String> sampleIds = vcfReader.vcfHeader().getGenotypeSamples();
        int tumorOrdinal = -1;

        List<String> genotypeSamples = vcfReader.vcfHeader().getGenotypeSamples();

        for(int genotypeIndex = 0; genotypeIndex < genotypeSamples.size(); ++genotypeIndex)
        {
            String genotypeSample = genotypeSamples.get(genotypeIndex);

            if(genotypeSample.equals(sampleId))
            {
                tumorOrdinal = genotypeIndex;
                break;
            }
        }

        if(tumorOrdinal == -1)
        {
            LNX_LOGGER.error("vcf({}) missing sampleId(sampleId)", vcfFile, sampleId);
            return Collections.emptyList();
        }

        int referenceOrdinal = -1;

        if(sampleIds.size() > 1)
        {
            referenceOrdinal = tumorOrdinal == 0 ? 1 : 0;
        }

        svFactory.setGenotypeOrdinals(referenceOrdinal, tumorOrdinal);

        for(VariantContext variantContext : vcfReader.iterator())
        {
            svFactory.addVariantContext(variantContext);
        }

        List<StructuralVariantData> svDataList = Lists.newArrayList();
        int svId = 0;

        for(StructuralVariant var : svFactory.results())
        {
            svDataList.add(convertGermlineSvData(var, svId++));
        }

        LNX_LOGGER.info("loaded {} germline SV data records from VCF file: {}", svDataList.size(), vcfFile);
        return svDataList;
    }

    public static List<SvVarData> createSvData(final List<StructuralVariantData> svRecords, final LinxConfig config)
    {
        List<SvVarData> svDataItems = Lists.newArrayList();

        for(final StructuralVariantData svRecord : svRecords)
        {
            final String filter = svRecord.filter();

            if(filter.isEmpty() || filter.equals(PASS) || filter.equals(INFERRED) || config.IsGermline)
            {
                svDataItems.add(new SvVarData(svRecord));
            }
        }

        return svDataItems;
    }

    private static double extractPurpleArrayValue(final VariantContext context, final String vcfTag)
    {
        if(context == null)
            return 0;

        List<Double> array = context.getAttributeAsDoubleList(vcfTag, 0);
        return array != null && !array.isEmpty() ? array.get(0) : 0;
    }

    public static StructuralVariantData convertGermlineSvData(final StructuralVariant var, int svId)
    {
        final VariantContext contextStart = var.startContext();
        final VariantContext contextEnd = var.endContext();

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
                .junctionCopyNumber(contextStart.getAttributeAsDouble(PURPLE_JUNCTION_COPY_NUMBER, 0))
                .startAF(valueNotNull(var.start().alleleFrequency()))
                .endAF(var.end() == null ? 0 : valueNotNull(var.end().alleleFrequency()))
                .adjustedStartAF(extractPurpleArrayValue(contextStart, PURPLE_AF))
                .adjustedEndAF(extractPurpleArrayValue(contextEnd, PURPLE_AF))
                .adjustedStartCopyNumber(extractPurpleArrayValue(contextStart, PURPLE_CN))
                .adjustedEndCopyNumber(extractPurpleArrayValue(contextEnd, PURPLE_CN))
                .adjustedStartCopyNumberChange(extractPurpleArrayValue(contextStart, PURPLE_CN_CHANGE))
                .adjustedEndCopyNumberChange(extractPurpleArrayValue(contextEnd, PURPLE_CN_CHANGE))
                .insertSequence(var.insertSequence())
                .type(var.type())
                .filter(var.filter())
                .imprecise(var.imprecise())
                .qualityScore(valueNotNull(var.qualityScore()))
                .event(valueNotNull(var.event()))
                .startTumorVariantFragmentCount(valueNotNull(var.start().tumorVariantFragmentCount()))
                .startTumorReferenceFragmentCount(valueNotNull(var.start().tumorReferenceFragmentCount()))
                .startNormalVariantFragmentCount(valueNotNull(var.start().normalVariantFragmentCount()))
                .startNormalReferenceFragmentCount(valueNotNull(var.start().normalReferenceFragmentCount()))
                .endTumorVariantFragmentCount(var.end() == null ? 0 : valueNotNull(var.end().tumorVariantFragmentCount()))
                .endTumorReferenceFragmentCount(var.end() == null ? 0 : valueNotNull(var.end().tumorReferenceFragmentCount()))
                .endNormalVariantFragmentCount(var.end() == null ? 0 : valueNotNull(var.end().normalVariantFragmentCount()))
                .endNormalReferenceFragmentCount(var.end() == null ? 0 : valueNotNull(var.end().normalReferenceFragmentCount()))
                .startIntervalOffsetStart(valueNotNull(var.start().startOffset()))
                .startIntervalOffsetEnd(valueNotNull(var.start().endOffset()))
                .endIntervalOffsetStart(var.end() == null ? 0 : valueNotNull(var.end().startOffset()))
                .endIntervalOffsetEnd(var.end() == null ? 0 : valueNotNull(var.end().endOffset()))
                .inexactHomologyOffsetStart(valueNotNull(var.start().inexactHomologyOffsetStart()))
                .inexactHomologyOffsetEnd(valueNotNull(var.start().inexactHomologyOffsetEnd()))
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
                .insertSequenceRepeatOrientation(valueNotNull(var.insertSequenceRepeatOrientation()))
                .insertSequenceRepeatCoverage(valueNotNull(var.insertSequenceRepeatCoverage()))
                .startAnchoringSupportDistance(var.start().anchoringSupportDistance())
                .endAnchoringSupportDistance(var.end() == null ? 0 : var.end().anchoringSupportDistance())
                .ponCount(var.startContext().getAttributeAsInt(PON_COUNT, 0))
                .build();
    }

    private static double valueNotNull(final Double value)
    {
        return value != null ? value : 0D;
    }
    private static int valueNotNull(final Integer value)
    {
        return value != null ? value : 0;
    }
    private static byte valueNotNull(final Byte value)
    {
        return value != null ? value : 0;
    }
    private static String valueNotNull(final String value)
    {
        return value != null ? value : Strings.EMPTY;
    }

}
