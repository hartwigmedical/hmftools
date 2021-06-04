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
import com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleles;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.read.SAMRecordReader;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_CHR;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_GENES;

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
    private final List<Integer> mVariantLoci;
    private final Map<String,Map<Integer,Set<String>>> mHetLociSansVariants;
    private final LilacConfig mConfig;

    private final List<SomaticVariant> mSomaticVariants;

    private final Set<CodingEffect> UNKNOWN_CODING_EFFECT;
    private final Map<String,TranscriptData> mHlaTranscriptData;

    private final LociPosition mLociPositionFinder;

    public SomaticVariantAnnotation(
            final LilacConfig config, final Map<String, TranscriptData> transcriptData,
            final Map<String, Map<Integer, Set<String>>> geneAminoAcidHetLociMap, final LociPosition lociPositionFinder)
    {
        mConfig = config;
        mHlaTranscriptData = transcriptData;
        UNKNOWN_CODING_EFFECT = Sets.newHashSet(CodingEffect.NONE, CodingEffect.UNDEFINED);

        mHetLociSansVariants = Maps.newHashMap();
        mVariantLoci = Lists.newArrayList();

        mLociPositionFinder = lociPositionFinder;

        mSomaticVariants = Lists.newArrayList();
        loadSomaticVariants();

        if(!mSomaticVariants.isEmpty())
            buildLoci(geneAminoAcidHetLociMap);
    }

    public List<SomaticVariant> getSomaticVariants() { return mSomaticVariants; }

    private void buildLoci(final Map<String, Map<Integer, Set<String>>> geneAminoAcidHetLociMap)
    {
        mSomaticVariants.stream()
                .map(x -> mLociPositionFinder.nucelotideLoci(x.Position))
                .filter(x -> x >= 0)
                .mapToInt(x -> x / 3)
                .forEach(x -> mVariantLoci.add(x));

        for(Map.Entry<String,Map<Integer,Set<String>>> geneEntry : geneAminoAcidHetLociMap.entrySet())
        {
            Map<Integer,Set<String>> lociSeqMap = geneEntry.getValue().entrySet().stream()
                    .filter(x -> !mVariantLoci.contains(x.getKey()))
                    .collect(Collectors.toMap(entry -> entry.getKey(), entry -> entry.getValue()));

            mHetLociSansVariants.put(geneEntry.getKey(), lociSeqMap);
        }
    }

    public final List<HlaAlleleCoverage> assignAlleleCoverage(
            final SomaticVariant variant, final SAMRecordReader reader, final List<HlaSequenceLoci> winners)
    {
        List<Fragment> fragments = reader.readFromBam(variant);

        fragments.forEach(x -> x.qualityFilter(mConfig.MinBaseQual));
        fragments.forEach(x -> x.buildAminoAcids());

        FragmentAlleleMapper fragAlleleMapper = new FragmentAlleleMapper(mHetLociSansVariants, Maps.newHashMap(), Lists.newArrayList());

        List<FragmentAlleles> variantFragmentAlleles = fragAlleleMapper.createFragmentAlleles(fragments, winners, Lists.newArrayList());

        List<HlaAlleleCoverage> coverage = HlaAlleleCoverage.proteinCoverage(variantFragmentAlleles);

        Collections.sort(coverage, new HlaAlleleCoverage.TotalCoverageSorter());

        return coverage.stream().filter(x -> x.TotalCoverage == coverage.get(0).TotalCoverage).collect(Collectors.toList());
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
