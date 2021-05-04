package com.hartwig.hmftools.lilac.out;

import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;
import com.hartwig.hmftools.lilackt.cna.HlaCopyNumber;
import com.hartwig.hmftools.lilackt.hla.HlaAllele;

import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;

import org.jetbrains.annotations.NotNull;

public class HlaOut
{
    private final HlaComplexCoverage mReferenceCoverage;
    private final HlaComplexCoverage tumorCoverage;
    private final List<HlaCopyNumber> tumorCopyNumber;
    private final List<SomaticCodingCount> somaticCodingCount;

    public HlaOut(final HlaComplexCoverage referenceCoverage, final HlaComplexCoverage tumorCoverage,
            final List<HlaCopyNumber> tumorCopyNumber, final List<SomaticCodingCount> somaticCodingCount)
    {
        this.mReferenceCoverage = referenceCoverage;
        this.tumorCoverage = tumorCoverage;
        this.tumorCopyNumber = tumorCopyNumber;
        this.somaticCodingCount = somaticCodingCount;
    }

    public final void write(final String fileName)
    {
        /*
        File file = new File(fileName);
        FilesKt.writeText$default((File) file, (String) (this.generateAlleleHeader() + "\n"), null, (int) 2, null);
        int n = 0;
        int n2 = 5;
        while(n <= n2)
        {
            void i;
            FilesKt.appendText$default((File) file, (String) (this.generateAlleleBody((int) i) + "\n"), null, (int) 2, null);
            ++i;
        }
        
         */
    }

    private final String generateAlleleHeader()
    {
        StringJoiner header = new StringJoiner("\t").add("Allele")
                .add("RefTotal")
                .add("RefUnique")
                .add("RefShared")
                .add("RefWild")
                .add("TumorTotal")
                .add("TumorUnique")
                .add("TumorShared")
                .add("TumorWild")
                .add("TumorCopyNumber")
                .add("SomaticMissense")
                .add("SomaticNonsenseOrFrameshift")
                .add("SomaticSplice")
                .add("SomaticSynonymous")
                .add("SomaticInframeIndel");
        String string = header.toString();
        return string;
    }

    private final String generateAlleleBody(int index)
    {
        /*
        generateAlleleBody .1 $fun$format$1 = generateAlleleBody .1.INSTANCE;
        HlaAlleleCoverage ref = this.referenceCoverage.getAlleleCoverage().get(index);
        Collection collection = this.tumorCoverage.getAlleleCoverage();
        HlaAlleleCoverage tumor = !collection.isEmpty()
                ? this.tumorCoverage.getAlleleCoverage().get(index)
                : new HlaAlleleCoverage(ref.getAllele(), 0, 0.0, 0.0);
        double copyNumber = this.tumorCopyNumber.get(index).getCopyNumber();
        SomaticCodingCount codingCount = this.somaticCodingCount.get(index);
        double d = ref.getTotalCoverage();
        StringJoiner stringJoiner = new StringJoiner("\t").add(ref.getAllele().asFourDigit().toString())
                .add(String.valueOf((int) Math.rint(d)))
                .add(String.valueOf(ref.getUniqueCoverage()));
        d = ref.getSharedCoverage();
        StringJoiner stringJoiner2 = stringJoiner.add(String.valueOf((int) Math.rint(d)));
        d = ref.getWildCoverage();
        StringJoiner stringJoiner3 = stringJoiner2.add(String.valueOf((int) Math.rint(d)));
        d = tumor.getTotalCoverage();
        StringJoiner stringJoiner4 = stringJoiner3.add(String.valueOf((int) Math.rint(d))).add(String.valueOf(tumor.getUniqueCoverage()));
        d = tumor.getSharedCoverage();
        StringJoiner stringJoiner5 = stringJoiner4.add(String.valueOf((int) Math.rint(d)));
        d = tumor.getWildCoverage();
        StringJoiner header = stringJoiner5.add(String.valueOf((int) Math.rint(d)))
                .add($fun$format$1.invoke(copyNumber, 2))
                .add($fun$format$1.invoke(codingCount.getMissense(), 2))
                .add($fun$format$1.invoke(codingCount.getNonsense(), 2))
                .add($fun$format$1.invoke(codingCount.getSplice(), 2))
                .add($fun$format$1.invoke(codingCount.getSynonymous(), 2))
                .add($fun$format$1.invoke(codingCount.getInframeIndel(), 2));
        String string = header.toString();
        Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "header.toString()");
        return string;
        
         */
        
        return "";
    }

    public static HlaOut create(final HlaComplexCoverage referenceCoverage, final HlaComplexCoverage tumorCoverage,
            final List<HlaCopyNumber> tumorCopyNumber, final List<SomaticCodingCount> somaticCodingCount)
    {
        return null;

        /*
        boolean bl;
        Iterable $receiver$iv;
        Intrinsics.checkParameterIsNotNull((Object) referenceCoverage, (String) "referenceCoverage");
        Intrinsics.checkParameterIsNotNull((Object) tumorCoverage, (String) "tumorCoverage");
        Intrinsics.checkParameterIsNotNull(tumorCopyNumber, (String) "tumorCopyNumber");
        Intrinsics.checkParameterIsNotNull(somaticCodingCount, (String) "somaticCodingCount");
        Iterable iterable = $receiver$iv = (Iterable) tumorCopyNumber;
        Object object = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                HlaCopyNumber it = (HlaCopyNumber) a;
                boolean bl = false;
                Comparable comparable = it.getAllele();
                it = (HlaCopyNumber) b;
                Comparable comparable2 = comparable;
                bl = false;
                HlaAllele hlaAllele = it.getAllele();
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) hlaAllele);
            }
        };
        List sortedCopyNumber = CollectionsKt.sortedWith((Iterable) iterable, (Comparator) object);
        Iterable $receiver$iv2 = somaticCodingCount;
        object = $receiver$iv2;
        Comparator comparator = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                SomaticCodingCount it = (SomaticCodingCount) a;
                boolean bl = false;
                Comparable comparable = it.getAllele();
                it = (SomaticCodingCount) b;
                Comparable comparable2 = comparable;
                bl = false;
                HlaAllele hlaAllele = it.getAllele();
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) hlaAllele);
            }
        };
        List sortedCodingCount = CollectionsKt.sortedWith((Iterable) object, (Comparator) comparator);
        boolean bl2 = bl = sortedCopyNumber.size() == somaticCodingCount.size();
        if(!bl)
        {
            object = "Failed requirement.";
            throw (Throwable) new IllegalArgumentException(object.toString());
        }
        boolean bl3 = bl = referenceCoverage.getAlleleCoverage().size() == sortedCopyNumber.size();
        if(!bl)
        {
            object = "Failed requirement.";
            throw (Throwable) new IllegalArgumentException(object.toString());
        }
        return new HlaOut(referenceCoverage, tumorCoverage, sortedCopyNumber, sortedCodingCount);

         */
    }
}
