package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.gripss.GripssConfig.GR_LOGGER;

public class VcfWriter
{
    private final GripssConfig mConfig;

    public VcfWriter(final GripssConfig config)
    {
        mConfig = config;
    }

    public void write()
    {
        String outputFile = mConfig.OutputDir + "gripss.filtered.vcf.gz";
        GR_LOGGER.info("writing output VCF file: {}", outputFile);

        /*
        logger.info("Writing file: ${config.outputVcf}")
        val combinedLinks = LinkStore(combinedTransitiveAssemblyLinks, dsbLinks)
        val finalFilters: SoftFilterStore = softFiltersAfterSingleDedup.update(setOf(), allRescues)
        fileWriter.writeHeader(version.version(), fileReader.fileHeader, outputSampleNames)
        for (variant in variantStore.selectAll()) {

            val localLinkedBy = combinedLinks[variant.vcfId]
            val remoteLinkedBy = combinedLinks[variant.mateId]
            val altPath = alternatePathsStringsByVcfId[variant.vcfId]

            val filters = finalFilters.filters(variant.vcfId, variant.mateId)
            fileWriter.writeVariant(variant.context(localLinkedBy, remoteLinkedBy, altPath, hotspots.contains(variant.vcfId), filters))
        }
        */

        // combine transitive, assembly and DSB links into a single set - where combining is done by VCF id, or link pair IDs??



        // write variants in chromosome & position order, ie using chromosome map and breakend  list

    }

    /*
        fun context(localLink: String, remoteLink: String, altPath: String?, isHotspot: Boolean, filters: Set<String>): VariantContext {
        val genotypesToWrite = mutableListOf(tumorGenotype)
        normalGenotype?.let { x -> genotypesToWrite.add(x) }

        // genotypes are of type FastGenotype, which is a 3rd-party class

        // need to record if the SV matches a known fusion pair

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

}
