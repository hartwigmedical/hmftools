package com.hartwig.hmftools.lilac.variant;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_GENES;
import static com.hartwig.hmftools.lilac.LilacConstants.shortGeneName;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.logging.log4j.core.util.FileUtils;

public class CopyNumberAssignment
{
    private final Map<String,List<CopyNumberData>> mSampleCopyNumberData;

    public CopyNumberAssignment()
    {
        mSampleCopyNumberData = Maps.newHashMap();
    }

    public void loadCopyNumberData(final LilacConfig config)
    {
        if(config.CopyNumberFile.isEmpty())
            return;

        try
        {
            if(config.CopyNumberFile.contains(config.Sample))
            {
                List<GeneCopyNumber> hlaGeneCopyNumbers = GeneCopyNumberFile.read(config.CopyNumberFile).stream()
                        .filter(x -> HLA_GENES.contains(x.gene())).collect(Collectors.toList());

                List<CopyNumberData> cnDataList = hlaGeneCopyNumbers.stream()
                        .map(x -> new CopyNumberData(x.gene(), x.minCopyNumber(), x.minMinorAlleleCopyNumber()))
                        .collect(Collectors.toList());

                mSampleCopyNumberData.put(config.Sample, cnDataList);
            }
            else
            {
                // load a cohort file - for now only retain the required sample's data
                final List<String> fileData = Files.readAllLines(new File(config.CopyNumberFile).toPath());
                String header = fileData.get(0);
                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);
                fileData.remove(0); // remove header
                int sampleIndex = fieldsIndexMap.get("SampleId");
                int geneIndex = fieldsIndexMap.get("Gene");
                int minorIndex = fieldsIndexMap.get("MinCopyNumber");
                int minorMinIndex = fieldsIndexMap.get("MinMinorAlleleCopyNumber");

                boolean loadAll = false;

                if(loadAll)
                {
                    List<CopyNumberData> cnDataList = null;
                    String currentSample = "";

                    for(final String line : fileData)
                    {
                        String[] items = line.split(DELIM, -1);

                        String sampleId = items[sampleIndex];

                        if(!currentSample.equals(sampleId))
                        {
                            cnDataList = Lists.newArrayList();
                            currentSample = sampleId;
                            mSampleCopyNumberData.put(sampleId, cnDataList);
                        }

                        cnDataList.add(new CopyNumberData(
                                items[geneIndex], Double.parseDouble(items[minorIndex]), Double.parseDouble(items[minorMinIndex])));
                    }
                }
                else
                {
                    // only load the configured sample
                    List<CopyNumberData> cnDataList = Lists.newArrayList();

                    for(final String line : fileData)
                    {
                        String[] items = line.split(DELIM, -1);

                        String sampleId = items[sampleIndex];

                        if(!sampleId.equals(config.Sample))
                        {
                            if(!cnDataList.isEmpty())
                                break;
                            else
                                continue;
                        }

                        cnDataList.add(new CopyNumberData(
                                items[geneIndex], Double.parseDouble(items[minorIndex]), Double.parseDouble(items[minorMinIndex])));
                    }

                    mSampleCopyNumberData.put(config.Sample, cnDataList);
                }
            }
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to read gene copy number file({}): {}", config.CopyNumberFile, e.toString());
        }
    }

    public static Map<HlaAllele,Double> formEmptyAlleleCopyNumber(final List<HlaAllele> winners)
    {
        return winners.stream().collect(Collectors.toMap(entry -> entry, entry -> new Double(0)));
    }

    public Map<HlaAllele,Double> assign(
            final String sampleId, final List<HlaAllele> winners, final HlaComplexCoverage tumorCoverage)
    {
        LL_LOGGER.info("calculating tumor copy number of winning alleles");

        List<CopyNumberData> cnDataList = mSampleCopyNumberData.get(sampleId);

        if(cnDataList == null || cnDataList.isEmpty())
            return formEmptyAlleleCopyNumber(winners);

        Map<HlaAllele,Double> alleleCopyNumbers = Maps.newHashMap();

        for(String gene : HLA_GENES)
        {
            String geneId = shortGeneName(gene);

            CopyNumberData cnData = cnDataList.stream().filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

            if(cnData == null)
            {
                LL_LOGGER.warn("missing gene({}) copy number data", gene);
                continue;
            }

            List<HlaAlleleCoverage> coverage = tumorCoverage.getAlleleCoverage().stream()
                    .filter(x -> x.Allele.Gene.equals(geneId)).collect(Collectors.toList());

            if(coverage.size() != 2)
            {
                LL_LOGGER.warn("missing gene({}) allele coverage", gene);
                continue;
            }

            alleleCopyNumber(cnData, alleleCopyNumbers, coverage.get(0), coverage.get(1));
        }

        return alleleCopyNumbers;
    }

    private void alleleCopyNumber(
            final CopyNumberData copyNumberData, final Map<HlaAllele,Double> alleleCopyNumbers,
            final HlaAlleleCoverage alleleCoverage1, final HlaAlleleCoverage alleleCoverage2)
    {

        double minor = copyNumberData.MinMinorAlleleCopyNumber;
        double major = copyNumberData.MinCopyNumber - minor;

        if(alleleCoverage1.TotalCoverage >= alleleCoverage2.TotalCoverage)
        {
            alleleCopyNumbers.put(alleleCoverage1.Allele, major);
            alleleCopyNumbers.put(alleleCoverage2.Allele, minor);
        }
        else
        {
            alleleCopyNumbers.put(alleleCoverage2.Allele, major);
            alleleCopyNumbers.put(alleleCoverage1.Allele, minor);
        }
    }

}
