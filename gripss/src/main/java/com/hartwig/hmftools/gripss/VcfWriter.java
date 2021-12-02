package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_ALT_PATH;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_EVENT_TYPE;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_HOTSPOT;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_LOCAL_LINKED_BY;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REALIGN;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REMOTE_LINKED_BY;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_TAF;
import static com.hartwig.hmftools.gripss.filters.FilterType.HARD_FILTERED;
import static com.hartwig.hmftools.gripss.filters.FilterType.PON;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.filters.FilterType;
import com.hartwig.hmftools.gripss.links.LinkStore;

import htsjdk.variant.variantcontext.Genotype;
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
    private final GripssConfig mConfig;
    private final GenotypeIds mGenotypeIds;
    private final SvDataCache mDataCache;
    private final FilterCache mFilterCache;

    private final VariantContextWriter mUnfilteredWriter;
    private final VariantContextWriter mFilteredWriter;

    public VcfWriter(
            final GripssConfig config, final VCFHeader vcfHeader, final String gripssVersion, final GenotypeIds genotypeIds,
            final SvDataCache dataCache, final FilterCache filterCache)
    {
        mConfig = config;
        mGenotypeIds = genotypeIds;
        mFilterCache = filterCache;
        mDataCache = dataCache;

        final String suffix = config.OutputId != null ? "." + config.OutputId + ".vcf.gz" : ".vcf.gz";
        final String unfilteredVcf = config.OutputDir + config.SampleId + ".gripss" + suffix;
        final String filteredVcf = config.OutputDir + config.SampleId + ".gripss.filtered" + suffix;

        mFilteredWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(vcfHeader.getSequenceDictionary())
                .setOutputFile(filteredVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();

        mUnfilteredWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(vcfHeader.getSequenceDictionary())
                .setOutputFile(unfilteredVcf)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();

        writeHeader(mFilteredWriter, vcfHeader, gripssVersion);
        writeHeader(mUnfilteredWriter, vcfHeader, gripssVersion);
    }

    private void writeHeader(
            final VariantContextWriter writer, final VCFHeader vcfHeader, final String gripssVersion)
    {
        VCFHeader newHeader = new VCFHeader(vcfHeader);
        newHeader.addMetaDataLine(new VCFHeaderLine("gripssVersion", gripssVersion));

        for(FilterType filter : FilterType.values())
        {
            if(filter == HARD_FILTERED)
                continue;

            newHeader.addMetaDataLine(new VCFFilterHeaderLine(
                    FilterType.vcfName(filter), FilterType.vcfInfoString(filter)));
        }

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(VT_REALIGN, 1, VCFHeaderLineType.Flag, "Variant was realigned"));
        newHeader.addMetaDataLine(new VCFInfoHeaderLine(VT_EVENT_TYPE, 1, VCFHeaderLineType.String, "Structural variant type"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                VT_TAF, 1, VCFHeaderLineType.Float,"Tumor allelic frequency (fragment support / total support)"));
        newHeader.addMetaDataLine(new VCFInfoHeaderLine(VT_ALT_PATH, 1, VCFHeaderLineType.String, "Alternate path"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                VT_LOCAL_LINKED_BY, 1, VCFHeaderLineType.String, "Breakend linking information"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                VT_REMOTE_LINKED_BY, 1, VCFHeaderLineType.String, "Partner breakend linking information"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(VT_HOTSPOT, 1, VCFHeaderLineType.Flag, "Variant is a hotspot"));

        // writer.writeHeader(VCFHeader(metaData, samples))
        writer.writeHeader(newHeader);
    }

    public void write(final LinkStore combinedLinks, final Map<Breakend,String> idPathMap)
    {
        for(HumanChromosome humanChromosome : HumanChromosome.values())
        {
            String chromosome = humanChromosome.toString();

            if(mConfig.RefGenVersion == V38)
                chromosome = RefGenomeFunctions.enforceChrPrefix(chromosome);

            List<Breakend> breakendList = mDataCache.getBreakendMap().get(chromosome);

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

                String altPathStr = idPathMap.get(breakend);

                writeBreakend(breakend, localLinks, remoteLinks, altPathStr);
            }

        }
        /*
        for (variant in variantStore.selectAll()) {

            val localLinkedBy = combinedLinks[variant.vcfId]
            val remoteLinkedBy = combinedLinks[variant.mateId]
            val altPath = alternatePathsStringsByVcfId[variant.vcfId]

            val filters = finalFilters.filters(variant.vcfId, variant.mateId)
            fileWriter.writeVariant(variant.context(localLinkedBy, remoteLinkedBy, altPath, hotspots.contains(variant.vcfId), filters))
        }
        */
    }

    private void writeBreakend(final Breakend breakend, final String localLinks, final String remoteLinks, final String altPathStr)
    {
        List<Genotype> genotypes = Lists.newArrayList(breakend.Context.getGenotype(mGenotypeIds.TumorOrdinal));

        if(mGenotypeIds.hasReference())
            genotypes.add(breakend.Context.getGenotype(mGenotypeIds.ReferenceOrdinal));

        VariantContextBuilder builder = new VariantContextBuilder(breakend.Context).genotypes(genotypes).filters();

        builder.log10PError(breakend.Qual / -10.0)
                .attribute(VT_TAF, String.format("%.4f", breakend.allelicFrequency()))
                .attribute(VT_HOTSPOT, mFilterCache.isHotspot(breakend.sv()))
                .attribute(VT_EVENT_TYPE, breakend.type());

        if(!localLinks.isEmpty())
            builder.attribute(VT_REMOTE_LINKED_BY, localLinks);

        if(!remoteLinks.isEmpty())
            builder.attribute(VT_LOCAL_LINKED_BY, remoteLinks);

        if(altPathStr != null && !altPathStr.isEmpty())
            builder.attribute(VT_ALT_PATH, altPathStr);

        List<FilterType> filters = mFilterCache.getBreakendFilters(breakend);

        if(filters == null)
            builder.filter(PASS);
        else
            filters.forEach(x -> builder.filter(FilterType.vcfName(x)));

        VariantContext variantContext = builder.make();

        mUnfilteredWriter.add(variantContext);

        if(filters.size() == 1 && (filters.contains(PASS) || filters.contains(PON)))
            mFilteredWriter.add(variantContext);
    }

    /*
        fun context(localLink: String, remoteLink: String, altPath: String?, isHotspot: Boolean, filters: Set<String>): VariantContext {
        val genotypesToWrite = mutableListOf(tumorGenotype)
        normalGenotype?.let { x -> genotypesToWrite.add(x) }

        val builder = VariantContextBuilder(context).genotypes(genotypesToWrite).filters()
        builder.log10PError(tumorQual / -10.0)
                .attribute(TAF, tumorAF)
                .attribute(LOCAL_LINKED_BY, localLink)
                .attribute(REMOTE_LINKED_BY, remoteLink)
                .attribute(HOTSPOT, isHotspot)
                .attribute(EVENTTYPE,variantType.eventType)

        altPath?.let { x -> builder.attribute(ALT_PATH, x) }
        filters.forEach { x -> builder.filter(x) }
        if (filters.isEmpty()) {
            builder.filter(PASS)
        }

        return builder.make()
    }

     */

    public void close()
    {
        mFilteredWriter.close();
        mUnfilteredWriter.close();
    }

}
