package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOTSPOT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOTSPOT_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.PON_COUNT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_COVERAGE_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_CLASS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_COVERAGE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_CLASS_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_TYPE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_TYPE_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.VCF_ZIP_EXTENSION;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.esvee.common.FileCommon.ESVEE_FILE_ID;
import static com.hartwig.hmftools.esvee.common.FileCommon.FILE_NAME_DELIM;
import static com.hartwig.hmftools.esvee.common.FileCommon.formOutputFile;
import static com.hartwig.hmftools.esvee.common.FilterType.PON;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.esvee.caller.annotation.RepeatMaskAnnotation;
import com.hartwig.hmftools.esvee.common.FilterType;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfWriter
{
    private final CallerConfig mConfig;
    private final GenotypeIds mGenotypeIds;
    private final SvDataCache mDataCache;

    private final VariantContextWriter mUnfilteredWriter;
    private final VariantContextWriter mSomaticWriter;
    private final VariantContextWriter mGermlineWriter;

    public VcfWriter(
            final CallerConfig config, final VCFHeader vcfHeader, final String gripssVersion, final GenotypeIds genotypeIds,
            final SvDataCache dataCache)
    {
        mConfig = config;
        mGenotypeIds = genotypeIds;
        mDataCache = dataCache;

        String fileSampleId = config.fileSampleId();

        String unfilteredVcf = formVcfFilename(fileSampleId, "unfiltered");
        mUnfilteredWriter = initialiseWriter(vcfHeader, gripssVersion, unfilteredVcf);

        if(mConfig.hasTumor())
        {
            String somaticVcf = formVcfFilename(fileSampleId, "somatic");
            mSomaticWriter = initialiseWriter(vcfHeader, gripssVersion, somaticVcf);
        }
        else
        {
            mSomaticWriter = null;
        }

        if(mConfig.hasReference())
        {
            String germlineVcf = formVcfFilename(fileSampleId, "germline");
            mGermlineWriter = initialiseWriter(vcfHeader, gripssVersion, germlineVcf);
        }
        else
        {
            mGermlineWriter = null;
        }
    }

    private String formVcfFilename(final String sampleId, final String fileId)
    {
        // write in format: SAMPLE_ID.esvee.run_id.somatic.vcf.gz (or without output ID if not in config)
        String outputId = mConfig.OutputId != null ? mConfig.OutputId + FILE_NAME_DELIM + fileId : fileId;
        return formOutputFile(mConfig.OutputDir, sampleId, ESVEE_FILE_ID, VCF_ZIP_EXTENSION.substring(1), outputId);
    }

    private VariantContextWriter initialiseWriter(final VCFHeader vcfHeader, final String esveeVersion, final String vcfFilename)
    {
        VariantContextWriter writer = new VariantContextWriterBuilder()
                .setReferenceDictionary(vcfHeader.getSequenceDictionary())
                .setOutputFile(vcfFilename)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();

        // ensure genotype sample IDs match the config - also done per variant
        List<String> genotypeSampleNames = Lists.newArrayList();

        if(mConfig.hasReference())
            genotypeSampleNames.add(mConfig.ReferenceId);

        if(mConfig.hasTumor())
            genotypeSampleNames.add(mConfig.TumorId);

        VCFHeader newHeader = new VCFHeader(vcfHeader.getMetaDataInInputOrder(), genotypeSampleNames);

        newHeader.addMetaDataLine(new VCFHeaderLine("esveeVersion", esveeVersion));

        newHeader.addMetaDataLine(new VCFFilterHeaderLine(PASS, "Variant passes all filters"));

        for(FilterType filter : FilterType.values())
        {
            newHeader.addMetaDataLine(new VCFFilterHeaderLine(filter.vcfTag(), filter.vcfDesc()));
        }

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(ALLELE_FRACTION, 1, VCFHeaderLineType.Float, ALLELE_FRACTION_DESC));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT, 1, VCFHeaderLineType.Flag, HOTSPOT_DESC));
        newHeader.addMetaDataLine(new VCFInfoHeaderLine(PON_COUNT, 1, VCFHeaderLineType.Integer, "PON count if in PON"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                REPEAT_MASK_REPEAT_CLASS, 1, VCFHeaderLineType.String, REPEAT_MASK_REPEAT_CLASS_DESC));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                REPEAT_MASK_REPEAT_TYPE, 1, VCFHeaderLineType.String, REPEAT_MASK_REPEAT_TYPE_DESC));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                REPEAT_MASK_COVERAGE, 1, VCFHeaderLineType.Float, REPEAT_MASK_COVERAGE_DESC));

        writer.writeHeader(newHeader);

        return writer;
    }

    public void writeBreakends()
    {
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            List<Breakend> breakendList = mDataCache.getBreakendMap().get(chrStr);

            if(breakendList == null)
                continue;

            for(Breakend breakend : breakendList)
            {
                writeBreakend(breakend);
            }
        }
    }

    private void writeBreakend(final Breakend breakend)
    {
        final Variant var = breakend.sv();

        List<Genotype> genotypes = Lists.newArrayList();

        if(mGenotypeIds.hasTumor())
            buildGenotype(genotypes, breakend, mConfig.TumorId);

        if(mGenotypeIds.hasReference())
            buildGenotype(genotypes, breakend, mConfig.ReferenceId);

        VariantContextBuilder builder = new VariantContextBuilder(breakend.Context).genotypes(genotypes).filters();

        // Qual = getGenotypeAttributeAsDouble(tumorGenotype, QUAL, 0);
        double qual = breakend.Context.getPhredScaledQual();
        builder.log10PError(qual / -10);

        Genotype tumorGenotype = genotypes.get(0);

        builder.attribute(ALLELE_FRACTION, breakend.calcAllelicFrequency(tumorGenotype));

        if(var.isHotspot())
            builder.attribute(HOTSPOT, true);

        if(var.ponCount() > 0)
            builder.attribute(PON_COUNT, breakend.sv().ponCount());

        if(var.getRmAnnotation() != null)
        {
            final RepeatMaskAnnotation rmAnnotation = var.getRmAnnotation();

            // remove any previously set
            builder.attribute(REPEAT_MASK_REPEAT_CLASS, rmAnnotation.RmData.ClassType);
            builder.attribute(REPEAT_MASK_REPEAT_TYPE, rmAnnotation.RmData.Repeat);
            builder.attribute(REPEAT_MASK_COVERAGE, rmAnnotation.Coverage);
        }

        // first write the unfiltered VCF with all breakends
        Set<FilterType> allFilters = breakend.sv().filters();

        writeBreakend(mUnfilteredWriter, builder, allFilters);

        if(breakend.sv().isGermline())
        {
            Set<FilterType> germlineFilters = allFilters.stream().filter(x -> !x.germlineOnly()).collect(Collectors.toSet());

            if(mGermlineWriter != null && germlineFilters.isEmpty())
                writeBreakend(mGermlineWriter, builder, germlineFilters);
        }
        else
        {
            Set<FilterType> somaticFilters = Sets.newHashSet(allFilters);

            if(mSomaticWriter != null)
            {
                if(somaticFilters.isEmpty() || (somaticFilters.size() == 1 && somaticFilters.contains(PON)))
                    writeBreakend(mSomaticWriter, builder, somaticFilters);
            }
        }
    }

    private void buildGenotype(final List<Genotype> genotypes, final Breakend breakend, final String sampleId)
    {
        Genotype existingGenotype = breakend.Context.getGenotype(sampleId);

        GenotypeBuilder genotypeBuilder = new GenotypeBuilder().copy(existingGenotype);

        genotypeBuilder.name(sampleId);

        int refPairSupport = getGenotypeAttributeAsInt(existingGenotype, REF_DEPTH_PAIR, 0);
        int refSupport = getGenotypeAttributeAsInt(existingGenotype, REF_DEPTH, 0);
        int varSupport = getGenotypeAttributeAsInt(existingGenotype, TOTAL_FRAGS, 0);
        int totalDepth = refPairSupport + refSupport + varSupport;

        genotypeBuilder.AD(new int[] { refPairSupport + refSupport, varSupport });
        genotypeBuilder.DP(totalDepth);

        genotypes.add(genotypeBuilder.make());
    }

    private static void writeBreakend(final VariantContextWriter writer, final VariantContextBuilder builder, final Set<FilterType> filters)
    {
        builder.getFilters().clear();

        if(filters.isEmpty())
            builder.filter(PASS);
        else
            filters.forEach(x -> builder.filter(x.vcfTag()));

        VariantContext variantContext = builder.make(true);

        writer.add(variantContext);
    }

    public void close()
    {
        mUnfilteredWriter.close();

        if(mSomaticWriter != null)
            mSomaticWriter.close();

        if(mGermlineWriter != null)
            mGermlineWriter.close();
    }
}
