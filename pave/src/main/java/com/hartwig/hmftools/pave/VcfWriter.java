package com.hartwig.hmftools.pave;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.VersionInfo.fromAppName;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact.VAR_TRANS_IMPACT_DELIM;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_COUNT;
import static com.hartwig.hmftools.common.variant.pon.PonCache.PON_MAX;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.APP_NAME;
import static com.hartwig.hmftools.pave.impact.PaveUtils.codonForBase;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.config.VersionInfo;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
import com.hartwig.hmftools.common.variant.pon.GnomadCache;
import com.hartwig.hmftools.common.variant.pon.PonCache;
import com.hartwig.hmftools.pave.annotation.Blacklistings;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;
import com.hartwig.hmftools.pave.annotation.Mappability;
import com.hartwig.hmftools.pave.annotation.ReferenceData;
import com.hartwig.hmftools.pave.annotation.Reportability;
import com.hartwig.hmftools.pave.impact.ProteinContext;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class VcfWriter
{
    private final VCFFileReader mHeader;
    private final VariantContextWriter mWriter;

    private final Map<HumanChromosome,List<VariantContext>> mChrPendingVariants;
    private final Set<HumanChromosome> mCompleteChromosomes;
    private HumanChromosome mCurrentChromosome;

    public boolean mWriteDetailed;

    protected static final String PROTEIN_CONTEXT = "PR_CXT";
    protected static final String PROTEIN_CONTEXT_DESC = "Protein context information";

    public VcfWriter(final String outputVCF, final String templateVCF)
    {
        mHeader = new VCFFileReader(new File(templateVCF), false);

        mWriter = new VariantContextWriterBuilder()
                .setReferenceDictionary(mHeader.getFileHeader().getSequenceDictionary())
                .setOutputFile(outputVCF)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                .build();

        mWriteDetailed = false;

        mChrPendingVariants = Maps.newHashMap();
        mCompleteChromosomes = Sets.newHashSet();
        mCurrentChromosome = HumanChromosome._1;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            mChrPendingVariants.put(chromosome, Lists.newArrayList());
        }
    }

    public void writeHeader(final ReferenceData referenceData, boolean setReportability, boolean writeDetailed)
    {
        final VersionInfo version = fromAppName(APP_NAME);

        VCFHeader newHeader = new VCFHeader(mHeader.getFileHeader());
        newHeader.addMetaDataLine(new VCFHeaderLine("PaveVersion", version.version()));

        VariantTranscriptImpact.writeHeader(newHeader);
        VariantImpactSerialiser.writeHeader(newHeader);

        if(writeDetailed)
        {
            mWriteDetailed = true;

            newHeader.addMetaDataLine(new VCFInfoHeaderLine(PROTEIN_CONTEXT, 1, VCFHeaderLineType.String, PROTEIN_CONTEXT_DESC));
        }

        if(referenceData.StandardPon.enabled() || referenceData.ArtefactsPon.enabled())
        {
            PonCache.addAnnotationHeader(newHeader);
        }

        PonCache.addFilterHeader(newHeader);

        if(referenceData.Gnomad.enabled())
        {
            GnomadCache.addAnnotationHeader(newHeader);
        }

        GnomadCache.addFilterHeader(newHeader);

        if(referenceData.VariantMappability.enabled())
        {
            Mappability.addHeader(newHeader);
        }

        if(referenceData.Clinvar.enabled())
        {
            ClinvarAnnotation.addHeader(newHeader);
        }

        if(referenceData.BlacklistedVariants.enabled())
        {
            Blacklistings.addHeader(newHeader);
        }

        if(setReportability)
            Reportability.addHeader(newHeader);

        mWriter.writeHeader(newHeader);
    }

    public VariantContext buildVariant(final VariantContext context, final VariantData variant, final VariantImpact variantImpact)
    {
        VariantContextBuilder builder = new VariantContextBuilder(variant.context())
                .genotypes(variant.context().getGenotypes())
                .filters(context.getFilters());

        if(!variant.filters().isEmpty())
        {
            builder.getFilters().addAll(variant.filters());
            builder.getFilters().remove(PASS_FILTER);
        }

        if(builder.getFilters().isEmpty())
            builder.getFilters().add(PASS_FILTER);

        VariantContext newContext = builder.make();

        if(!variant.getImpacts().isEmpty())
        {
            List<VariantTranscriptImpact> transImpacts = Lists.newArrayList();

            for(Map.Entry<String,List<VariantTransImpact>> entry : variant.getImpacts().entrySet())
            {
                final String geneName = entry.getKey();
                final List<VariantTransImpact> geneImpacts = entry.getValue();

                for(VariantTransImpact transImpact : geneImpacts)
                {
                    int affectedExon = transImpact.codingContext().ExonRank;
                    int affectedCodon = codonForBase(transImpact.codingContext().CodingBase);
                    String refSeqId = transImpact.TransData.RefSeqId == null ? "" : transImpact.TransData.RefSeqId;

                    transImpacts.add(new VariantTranscriptImpact(
                            transImpact.TransData.GeneId, geneName, transImpact.TransData.TransName,
                            transImpact.effectsStr(), transImpact.inSpliceRegion(),
                            transImpact.hgvsCoding(), transImpact.hgvsProtein(), refSeqId, affectedExon, affectedCodon));

                    if(mWriteDetailed && transImpact.TransData.IsCanonical)
                    {
                        writeProteinInfoString(newContext, transImpact);
                    }
                }
            }

            VariantTranscriptImpact.writeVcfData(newContext, transImpacts);

            if(variantImpact != null)
                VariantImpactSerialiser.writeImpactDetails(newContext, variantImpact);
        }

        if(variant.gnomadFrequency() != null && !newContext.getCommonInfo().hasAttribute(GNOMAD_FREQ))
            newContext.getCommonInfo().putAttribute(GNOMAD_FREQ, variant.gnomadFrequency());

        if((variant.ponSampleCount() > 0 || variant.ponMaxReadCount() > 0) && !newContext.getCommonInfo().hasAttribute(PON_COUNT))
        {
            newContext.getCommonInfo().putAttribute(PON_COUNT, variant.ponSampleCount());
            newContext.getCommonInfo().putAttribute(PON_MAX, variant.ponMaxReadCount());
        }

        return newContext;
    }

    private static void writeProteinInfoString(final VariantContext variantContext, final VariantTransImpact transImpact)
    {
        if(!transImpact.hasProteinContext() || variantContext.getCommonInfo().hasAttribute(PROTEIN_CONTEXT))
            return;

        final ProteinContext pc = transImpact.proteinContext();
        StringJoiner sj = new StringJoiner(VAR_TRANS_IMPACT_DELIM);
        sj.add(pc.RefCodonBasesExtended);
        sj.add(pc.AltCodonBasesComplete);
        sj.add(pc.RefAminoAcids);
        sj.add(pc.AltAminoAcids);
        sj.add(pc.NetRefAminoAcids);
        sj.add(pc.NetAltAminoAcids);
        sj.add(String.valueOf(pc.CodonIndex));
        sj.add(format("%d_%d", pc.NetCodonIndexRange[0], pc.NetCodonIndexRange[1]));

        variantContext.getCommonInfo().putAttribute(PROTEIN_CONTEXT, sj.toString());
    }

    public synchronized void onChromosomeComplete(final HumanChromosome chromosome)
    {
        mCompleteChromosomes.add(chromosome);
        writePendingVariants();
    }

    public synchronized void writeVariant(final HumanChromosome chromosome, final VariantContext variantContext)
    {
        if(mCurrentChromosome == chromosome)
        {
            mWriter.add(variantContext);
            return;
        }

        List<VariantContext> pendingVariants = mChrPendingVariants.get(chromosome);

        if(pendingVariants == null)
        {
            pendingVariants = Lists.newArrayList();
            mChrPendingVariants.put(chromosome, pendingVariants);
        }

        pendingVariants.add(variantContext);
    }

    private void writePendingVariants()
    {
        boolean exitOnNext = false;

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            if(!mCompleteChromosomes.contains(chromosome))
            {
                // move to the next chromosome to allow writing of current or pending variants
                mCurrentChromosome = chromosome;
                exitOnNext = true;
            }

            List<VariantContext> pendingVariants = mChrPendingVariants.get(chromosome);

            if(pendingVariants != null)
            {
                if(!pendingVariants.isEmpty())
                {
                    PV_LOGGER.debug("chr({}) writing {} pending variants", chromosome, pendingVariants.size());

                    pendingVariants.forEach(x -> mWriter.add(x));
                    pendingVariants.clear();
                }

                mChrPendingVariants.remove(chromosome);
            }

            if(exitOnNext)
                break;
        }
    }

    public void close()
    {
        writePendingVariants();
        mWriter.close();
    }
}
