package com.hartwig.hmftools.lilac.variant;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.LociPosition;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.read.BamReader;
import com.hartwig.hmftools.lilac.read.BamRecordReader;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_CHR;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_GENES;
import static com.hartwig.hmftools.lilac.seq.HlaSequence.WILD_STR;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class SomaticVariantAnnotation
{
    private final Map<String,List<Integer>> mGeneVariantLoci;
    private final LilacConfig mConfig;

    private final List<SomaticVariant> mSomaticVariants;

    private final Set<CodingEffect> UNKNOWN_CODING_EFFECT;
    private final Map<String,TranscriptData> mHlaTranscriptData;

    private final LociPosition mLociPositionFinder;

    public SomaticVariantAnnotation(
            final LilacConfig config, final Map<String, TranscriptData> transcriptData, final LociPosition lociPositionFinder)
    {
        mConfig = config;
        mHlaTranscriptData = transcriptData;
        UNKNOWN_CODING_EFFECT = Sets.newHashSet(CodingEffect.NONE, CodingEffect.UNDEFINED);

        mGeneVariantLoci = Maps.newHashMap();

        mLociPositionFinder = lociPositionFinder;

        mSomaticVariants = Lists.newArrayList();
        loadSomaticVariants();

        if(!mSomaticVariants.isEmpty())
        {
            for(SomaticVariant variant : mSomaticVariants)
            {
                int variantNucleotideLoci = mLociPositionFinder.nucelotideLoci(variant.Position);

                if(variantNucleotideLoci < 0)
                    continue;

                int variantAminoAcidLoci = variantNucleotideLoci / 3;

                List<Integer> geneLoci = mGeneVariantLoci.get(variant.Gene);

                if(geneLoci == null)
                {
                    geneLoci = Lists.newArrayList();
                    mGeneVariantLoci.put(variant.Gene, geneLoci);
                }

                geneLoci.add(variantAminoAcidLoci);
            }
        }
    }

    public List<SomaticVariant> getSomaticVariants() { return mSomaticVariants; }

    public final List<HlaAlleleCoverage> assignAlleleCoverage(
            final SomaticVariant variant, final BamReader reader, final List<HlaSequenceLoci> winners)
    {
        List<Fragment> fragments = reader.readFromBam(variant);
        fragments.forEach(x -> x.qualityFilter(mConfig.MinBaseQual));
        fragments.forEach(x -> x.buildAminoAcids());

        List<HlaAlleleCoverage> coverages = Lists.newArrayList();

        for(HlaSequenceLoci sequenceLoci : winners)
        {
            //if(!variant.Gene.equals(longGeneName(sequenceLoci.Allele.Gene)))
            //    continue;

            // any variant in this gene will be skipped
            List<Integer> variantLoci = mGeneVariantLoci.get(variant.Gene);

            int supportCount = 0;

            for(Fragment fragment : fragments)
            {
                boolean matches = true;
                int matchCount = 0;

                for(int locus = fragment.minAminoAcidLocus(); locus <= fragment.maxAminoAcidLocus(); ++locus)
                {
                    if(locus >= sequenceLoci.length())
                        break;

                    if(variantLoci.contains(locus))
                        continue;

                    int index = fragment.getAminoAcidLoci().indexOf(locus);
                    String fragmentAA = "";

                    if(index >= 0)
                    {
                        fragmentAA = fragment.getAminoAcids().get(index);
                    }
                    else
                    {
                        fragmentAA = fragment.getLowQualAminoAcid(locus);

                        if(fragmentAA.isEmpty())
                            continue;
                    }

                    String sequence = sequenceLoci.sequence(locus);

                    if(!fragmentAA.equals(sequence) && !sequence.equals(WILD_STR))
                    {
                        matches = false;
                        break;
                    }

                    ++matchCount;

                }

                if(matches && matchCount > 0)
                {
                    //LL_LOGGER.debug("allele({}) supported by fragment({} {}) matchedAAs({})",
                    //        seq.Allele, fragment.id(), fragment.readInfo(), matchCount);
                    ++supportCount;

                }
            }

            if(supportCount > 0)
                coverages.add(new HlaAlleleCoverage(sequenceLoci.Allele, supportCount, 0, 0));
        }

        // take top allele and any matching
        Collections.sort(coverages, new HlaAlleleCoverage.TotalCoverageSorter());
        return coverages.stream().filter(x -> x.TotalCoverage == coverages.get(0).TotalCoverage).collect(Collectors.toList());
    }

    private void loadSomaticVariants()
    {
        if(mConfig == null || mConfig.SomaticVariantsFile.isEmpty())
            return;

        if(mConfig.SomaticVariantsFile.contains(mConfig.Sample))
        {
            LL_LOGGER.info("reading somatic vcf: ", mConfig.SomaticVariantsFile);

            int minPosition = mHlaTranscriptData.values().stream().mapToInt(x -> x.TransStart).min().orElse(0);
            int maxPosition = mHlaTranscriptData.values().stream().mapToInt(x -> x.TransEnd).max().orElse(0);

            VCFFileReader fileReader = new VCFFileReader(new File(mConfig.SomaticVariantsFile), false);

            final CloseableIterator<VariantContext> variantIter = fileReader.isQueryable() ?
                    fileReader.query(HLA_CHR, minPosition, maxPosition) : fileReader.iterator();

            while(variantIter.hasNext())
            {
                VariantContext variant = variantIter.next();
                VariantContextDecorator enriched = new VariantContextDecorator(variant);

                if(HLA_GENES.contains(enriched.gene()) && enriched.isPass()
                        && !UNKNOWN_CODING_EFFECT.contains(enriched.canonicalCodingEffect()))
                {
                    mSomaticVariants.add(new SomaticVariant(
                            enriched.gene(), enriched.chromosome(), (int)enriched.position(), enriched.ref(), enriched.alt(),
                            enriched.filter(), enriched.canonicalCodingEffect(), enriched.context()));
                }
            }

            fileReader.close();
        }
        else
        {
            try
            {
                // load a cohort file - for now only retain the required sample's data
                final List<String> fileData = Files.readAllLines(new File(mConfig.SomaticVariantsFile).toPath());
                String header = fileData.get(0);
                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);
                fileData.remove(0); // remove header
                int sampleIndex = fieldsIndexMap.get("SampleId");
                int geneIndex = fieldsIndexMap.get("Gene");
                int chrIndex = fieldsIndexMap.get("Chromosome");
                int posIndex = fieldsIndexMap.get("Position");
                int refIndex = fieldsIndexMap.get("Ref");
                int altIndex = fieldsIndexMap.get("Alt");
                int filterIndex = fieldsIndexMap.get("Filter");
                int effectIndex = fieldsIndexMap.get("CanonicalCodingEffect");

                for(final String line : fileData)
                {
                    String[] items = line.split(DELIM, -1);

                    String sampleId = items[sampleIndex];

                    if(!sampleId.equals(mConfig.Sample))
                        continue;

                    String gene = items[geneIndex];
                    CodingEffect codingEffect = CodingEffect.valueOf(items[effectIndex]);
                    String filter = items[filterIndex];

                    if(!HLA_GENES.contains(gene))
                        continue;

                    if(UNKNOWN_CODING_EFFECT.contains(codingEffect))
                        continue;

                    if(!filter.equals(SomaticVariantFactory.PASS_FILTER))
                        continue;

                    mSomaticVariants.add(new SomaticVariant(
                            gene, items[chrIndex], Integer.parseInt(items[posIndex]), items[refIndex], items[altIndex],
                            filter, codingEffect, null));
                }
            }
            catch(IOException e)
            {
                LL_LOGGER.error("failed to read somatic variant file({}): {}", mConfig.SomaticVariantsFile, e.toString());
            }
        }

        LL_LOGGER.info("  found {} HLA variants", mSomaticVariants.size());
    }

}
