package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOTSPOT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOTSPOT_DESC;
import static com.hartwig.hmftools.common.sv.SvVcfTags.PON_COUNT;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_CLASS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_COVERAGE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REPEAT_MASK_REPEAT_TYPE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.VCF_ZIP_EXTENSION;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.esvee.common.FilterType.HARD_FILTERED;
import static com.hartwig.hmftools.esvee.common.FilterType.PON;
import static com.hartwig.hmftools.esvee.common.SvConstants.ESVEE_FILE_ID;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.esvee.caller.annotation.RepeatMaskAnnotation;
import com.hartwig.hmftools.esvee.common.FilterType;

import htsjdk.samtools.SAMSequenceRecord;
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
    private final FilterCache mFilterCache;

    // private final VariantContextWriter mUnfilteredWriter;
    private final VariantContextWriter mFilteredWriter;

    public VcfWriter(
            final CallerConfig config, final VCFHeader vcfHeader, final String gripssVersion, final GenotypeIds genotypeIds,
            final SvDataCache dataCache, final FilterCache filterCache)
    {
        mConfig = config;
        mGenotypeIds = genotypeIds;
        mFilterCache = filterCache;
        mDataCache = dataCache;

        String fileSampleId = config.GermlineMode && !config.ReferenceId.isEmpty() ? config.ReferenceId : config.SampleId;
        final String suffix = config.OutputId != null ? "." + config.OutputId + VCF_ZIP_EXTENSION : VCF_ZIP_EXTENSION;

        // final String unfilteredVcf = config.OutputDir + fileSampleId + ".gripss" + suffix;
        final String filteredVcf = config.OutputDir + fileSampleId + ESVEE_FILE_ID + suffix;

        mFilteredWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(vcfHeader.getSequenceDictionary())
                .setOutputFile(filteredVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();

        /*
        mUnfilteredWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(vcfHeader.getSequenceDictionary())
                .setOutputFile(unfilteredVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();
        */

        writeHeader(mFilteredWriter, vcfHeader, gripssVersion);
        // writeHeader(mUnfilteredWriter, vcfHeader, gripssVersion);
    }

    private void writeHeader(final VariantContextWriter writer, final VCFHeader vcfHeader, final String gripssVersion)
    {
        // ensure genotype sample IDs match the config - also done per variant
        List<String> genotypeSampleNames = Lists.newArrayList();

        if(!mConfig.ReferenceId.isEmpty())
            genotypeSampleNames.add(mConfig.ReferenceId);

        genotypeSampleNames.add(mConfig.SampleId);

        VCFHeader newHeader = new VCFHeader(vcfHeader.getMetaDataInInputOrder(), genotypeSampleNames);

        newHeader.addMetaDataLine(new VCFHeaderLine("gripssVersion", gripssVersion));

        for(FilterType filter : FilterType.values())
        {
            if(filter == HARD_FILTERED)
                continue;

            newHeader.addMetaDataLine(new VCFFilterHeaderLine(filter.vcfTag(), filter.vcfDesc()));
        }

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(HOTSPOT, 1, VCFHeaderLineType.Flag, HOTSPOT_DESC));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(PON_COUNT, 1, VCFHeaderLineType.Integer, "PON count if in PON"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                REPEAT_MASK_REPEAT_CLASS, 1, VCFHeaderLineType.String, "Inserted sequence repeatmasker repeat class"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                REPEAT_MASK_REPEAT_TYPE, 1, VCFHeaderLineType.String, "Inserted sequence repeatmasker repeat type"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                REPEAT_MASK_REPEAT_COVERAGE, 1, VCFHeaderLineType.Float,
                "Portion of inserted sequence whose alignment overlaps the repeatmasker repeat"));

        writer.writeHeader(newHeader);
    }

    public void write(final LinkStore combinedLinks, final VCFHeader vcfHeader)
    {
        final Map<SvData,List<FilterType>> svFiltersMap = Maps.newHashMap(); // to avoid collating SV filters twice

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(chromosome.toString());

            List<Breakend> breakendList = mDataCache.getBreakendMap().get(chrStr);

            if(breakendList == null)
                continue;

            for(Breakend breakend : breakendList)
            {
                String localLinks = combinedLinks.getBreakendLinksStr(breakend);
                String remoteLinks = "";

                Breakend otherBreakend = breakend.otherBreakend();

                if(otherBreakend != null)
                {
                    remoteLinks = combinedLinks.getBreakendLinksStr(otherBreakend);
                }

                writeBreakend(breakend, svFiltersMap, localLinks, remoteLinks);
            }
        }
    }

    private void writeBreakend(
            final Breakend breakend, final Map<SvData,List<FilterType>> svFiltersMap, final String localLinks, final String remoteLinks)
    {
        final SvData sv = breakend.sv();

        List<Genotype> genotypes = Lists.newArrayList();

        genotypes.add(new GenotypeBuilder().copy(breakend.Context.getGenotype(mGenotypeIds.TumorOrdinal)).name(mConfig.SampleId).make());

        if(mGenotypeIds.hasReference())
        {
            genotypes.add(new GenotypeBuilder().copy(breakend.Context.getGenotype(mGenotypeIds.ReferenceOrdinal)).name(mConfig.ReferenceId).make());
        }

        VariantContextBuilder builder = new VariantContextBuilder(breakend.Context).genotypes(genotypes).filters();

        builder.log10PError(breakend.Qual / -10.0);

        builder.attribute(HOTSPOT, mFilterCache.isHotspot(breakend.sv()));

        if(sv.ponCount() > 0)
            builder.attribute(PON_COUNT, breakend.sv().ponCount());

        if(sv.getRmAnnotation() != null)
        {
            final RepeatMaskAnnotation rmAnnotation = sv.getRmAnnotation();

            // remove any previously set
            builder.rmAttribute(REPEAT_MASK_REPEAT_CLASS);
            builder.rmAttribute(REPEAT_MASK_REPEAT_TYPE);
            builder.rmAttribute(REPEAT_MASK_REPEAT_COVERAGE);
            builder.attribute(REPEAT_MASK_REPEAT_CLASS, rmAnnotation.RmData.ClassType);
            builder.attribute(REPEAT_MASK_REPEAT_TYPE, rmAnnotation.RmData.Repeat);
            builder.attribute(REPEAT_MASK_REPEAT_COVERAGE, rmAnnotation.Coverage);
        }

        List<FilterType> svFilters;

        if(sv.isSgl())
        {
            svFilters = mFilterCache.getBreakendFilters(breakend);
        }
        else if(breakend == sv.breakendStart())
        {
            svFilters = mFilterCache.combineSvFilters(sv);
            svFiltersMap.put(sv, svFilters); // cache to avoid a second collation
        }
        else
        {
            svFilters = svFiltersMap.get(sv);
        }

        if(svFilters == null)
        {
            builder.filter(PASS);
        }
        else
        {
            svFilters.forEach(x -> builder.filter(x.vcfTag()));
        }

        VariantContext variantContext = builder.make(true);

        if(svFilters == null || (svFilters.size() == 1 && svFilters.get(0) == PON))
        {
            mFilteredWriter.add(variantContext);
        }

        // mUnfilteredWriter.add(variantContext);
    }

    public void close()
    {
        mFilteredWriter.close();
        // mUnfilteredWriter.close();
    }

}