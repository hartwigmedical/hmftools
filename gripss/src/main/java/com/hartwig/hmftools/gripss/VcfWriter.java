package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_ALT_PATH;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_EVENT_TYPE;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_HOTSPOT;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_LOCAL_LINKED_BY;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REALIGN;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REMOTE_LINKED_BY;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_RESCUE_INFO;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_TAF;
import static com.hartwig.hmftools.gripss.filters.FilterType.HARD_FILTERED;
import static com.hartwig.hmftools.gripss.filters.FilterType.PON;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterType;
import com.hartwig.hmftools.gripss.links.LinkRescue;
import com.hartwig.hmftools.gripss.links.LinkStore;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
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
                VT_LOCAL_LINKED_BY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Breakend linking information"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                VT_REMOTE_LINKED_BY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Partner breakend linking information"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(VT_HOTSPOT, 1, VCFHeaderLineType.Flag, "Variant is a hotspot"));

        newHeader.addMetaDataLine(new VCFInfoHeaderLine(
                VT_RESCUE_INFO, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Partner breakend rescue"));

        writer.writeHeader(newHeader);
    }

    public void write(
            final LinkStore combinedLinks, final Map<Breakend,String> idPathMap, final VCFHeader vcfHeader, final LinkRescue linkRescue)
    {
        final Map<SvData,List<FilterType>> svFiltersMap = Maps.newHashMap(); // to avoid collating SV filters twice

        for(SAMSequenceRecord seqRecord : vcfHeader.getSequenceDictionary().getSequences())
        {
            String chromosome = seqRecord.getContig();

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
                String rescueInfo = linkRescue.getRescueInfo().get(breakend);

                writeBreakend(breakend, svFiltersMap, localLinks, remoteLinks, altPathStr, rescueInfo);
            }
        }
    }

    private void writeBreakend(
            final Breakend breakend, final Map<SvData,List<FilterType>> svFiltersMap,
            final String localLinks, final String remoteLinks, final String altPathStr, String rescueInfo)
    {
        List<Genotype> genotypes = Lists.newArrayList(breakend.Context.getGenotype(mGenotypeIds.TumorOrdinal));

        if(mGenotypeIds.hasReference())
            genotypes.add(breakend.Context.getGenotype(mGenotypeIds.ReferenceOrdinal));

        VariantContextBuilder builder = new VariantContextBuilder(breakend.Context).genotypes(genotypes).filters();

        builder.log10PError(breakend.Qual / -10.0)
                .attribute(VT_TAF, String.format("%.4f", breakend.allelicFrequency()))
                .attribute(VT_HOTSPOT, mFilterCache.isHotspot(breakend.sv()))
                .attribute(VT_EVENT_TYPE, breakend.type());

        builder.rmAttribute(VT_ALT_PATH); // remove if set from an earlier run's file

        if(!localLinks.isEmpty())
            builder.attribute(VT_LOCAL_LINKED_BY, localLinks);
        else
            builder.attribute(VT_LOCAL_LINKED_BY, "");

        if(!remoteLinks.isEmpty())
            builder.attribute(VT_REMOTE_LINKED_BY, remoteLinks);
        else
            builder.attribute(VT_REMOTE_LINKED_BY, "");

        final SvData sv = breakend.sv();

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
            svFilters.forEach(x -> builder.filter(FilterType.vcfName(x)));
        }

        VariantContext variantContext = builder.make(true);

        if(svFilters == null || (svFilters.size() == 1 && svFilters.get(0) == PON))
        {
            mFilteredWriter.add(variantContext);
        }

        // write additional status and working data to unfiltered VCF
        if(altPathStr != null && !altPathStr.isEmpty())
            builder.attribute(VT_ALT_PATH, altPathStr);

        if(rescueInfo != null)
            builder.attribute(VT_RESCUE_INFO, rescueInfo);

        mUnfilteredWriter.add(variantContext);
    }

    public void close()
    {
        mFilteredWriter.close();
        mUnfilteredWriter.close();
    }

}
