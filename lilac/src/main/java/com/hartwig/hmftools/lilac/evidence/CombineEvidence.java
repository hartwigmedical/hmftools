package com.hartwig.hmftools.lilac.evidence;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

public class CombineEvidence
{
    // TODO
    public static final CombineEvidence INSTANCE = null;

    public final PhasedEvidence combine(@NotNull PhasedEvidence left, @NotNull PhasedEvidence common, @NotNull PhasedEvidence right)
    {
        PhasedEvidence leftWithCommon = this.combineOverlapping(left, common);
        PhasedEvidence result = this.combineOverlapping(leftWithCommon, right);
        return result;
    }

    public final boolean canCombine(@NotNull PhasedEvidence left, @NotNull PhasedEvidence common, @NotNull PhasedEvidence right)
    {
        return this.canCombine(left, common) && this.canCombine(common, right);
    }

    public static boolean canCombine(@NotNull PhasedEvidence left, @NotNull PhasedEvidence right)
    {
        return true;

        /*
        void $receiver$iv$iv;
        String string;
        Collection collection;
        void $receiver$iv$iv2;
        int it;
        Object element$iv;
        Iterator $receiver$iv$iv3;
        void $receiver$iv$iv4;
        int[] $receiver$iv;
        Intrinsics.checkParameterIsNotNull((Object) left, (String) "left");
        Intrinsics.checkParameterIsNotNull((Object) right, (String) "right");
        List indexIntersection =
                CollectionsKt.sorted((Iterable) CollectionsKt.toList((Iterable) CollectionsKt.intersect((Iterable) ArraysKt.toSet((int[]) left
                        .getAminoAcidIndices()), (Iterable) ArraysKt.toSet((int[]) right.getAminoAcidIndices()))));
        if(indexIntersection.isEmpty())
        {
            return false;
        }
        int[] nArray = $receiver$iv = left.getAminoAcidIndices();
        Object destination$iv$iv = new ArrayList();
        for(int element$iv$iv : $receiver$iv$iv4)
        {
            void it2 = element$iv$iv;
            boolean bl = false;
            if(!(!indexIntersection.contains((int) it2)))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        List uniqueToLeft = (List) destination$iv$iv;
        Object $receiver$iv2 = right.getAminoAcidIndices();
        destination$iv$iv = $receiver$iv2;
        Collection destination$iv$iv2 = new ArrayList();
        int n = ((void) $receiver$iv$iv3).length;
        for(int element$iv$iv = 0; element$iv$iv < n; ++element$iv$iv)
        {
            void element$iv$iv2;
            void it3 = element$iv$iv2 = $receiver$iv$iv3[element$iv$iv];
            boolean bl = false;
            if(!(!indexIntersection.contains((int) it3)))
            {
                continue;
            }
            destination$iv$iv2.add((int) element$iv$iv2);
        }
        List uniqueToRight = (List) destination$iv$iv2;
        $receiver$iv2 = uniqueToRight;
        if(!$receiver$iv2.isEmpty())
        {
            boolean bl;
            block19:
            {
                $receiver$iv2 = uniqueToLeft;
                if($receiver$iv2 instanceof Collection && ((Collection) $receiver$iv2).isEmpty())
                {
                    bl = false;
                }
                else
                {
                    $receiver$iv$iv3 = $receiver$iv2.iterator();
                    while($receiver$iv$iv3.hasNext())
                    {
                        element$iv = $receiver$iv$iv3.next();
                        it = ((Number) element$iv).intValue();
                        boolean bl2 = false;
                        Comparable comparable = CollectionsKt.min((Iterable) uniqueToRight);
                        if(comparable == null)
                        {
                            Intrinsics.throwNpe();
                        }
                        if(!(it > ((Number) ((Object) comparable)).intValue()))
                        {
                            continue;
                        }
                        bl = true;
                        break block19;
                    }
                    bl = false;
                }
            }
            if(bl)
            {
                return false;
            }
        }
        if(!($receiver$iv2 = (Collection) uniqueToLeft).isEmpty())
        {
            boolean bl;
            block20:
            {
                $receiver$iv2 = uniqueToRight;
                if($receiver$iv2 instanceof Collection && ((Collection) $receiver$iv2).isEmpty())
                {
                    bl = false;
                }
                else
                {
                    $receiver$iv$iv3 = $receiver$iv2.iterator();
                    while($receiver$iv$iv3.hasNext())
                    {
                        element$iv = $receiver$iv$iv3.next();
                        it = ((Number) element$iv).intValue();
                        boolean bl3 = false;
                        Comparable comparable = CollectionsKt.max((Iterable) uniqueToLeft);
                        if(comparable == null)
                        {
                            Intrinsics.throwNpe();
                        }
                        if(!(it < ((Number) ((Object) comparable)).intValue()))
                        {
                            continue;
                        }
                        bl = true;
                        break block20;
                    }
                    bl = false;
                }
            }
            if(bl)
            {
                return false;
            }
        }
        Iterable $receiver$iv3 = left.getEvidence().keySet();
        element$iv = $receiver$iv3;
        Iterable destination$iv$iv3 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv3, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv2)
        {
            String it3 = (String) item$iv$iv;
            collection = destination$iv$iv3;
            boolean bl = false;
            String $i$f$filterTo = it3;
            int n2 = left.getAminoAcidIndices().length - indexIntersection.size();
            String string2 = $i$f$filterTo;
            if(string2 == null)
            {
                throw new TypeCastException("null cannot be cast to non-null type java.lang.String");
            }
            Intrinsics.checkExpressionValueIsNotNull((Object) string2.substring(n2), (String) "(this as java.lang.String).substring(startIndex)");
            collection.add(string);
        }
        Set leftCommon = CollectionsKt.toSet((Iterable) ((List) destination$iv$iv3));
        Iterable $receiver$iv4 = right.getEvidence().keySet();
        destination$iv$iv3 = $receiver$iv4;
        Collection destination$iv$iv4 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv4, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it4;
            String bl = (String) item$iv$iv;
            collection = destination$iv$iv4;
            boolean bl4 = false;
            void var15_25 = it4;
            int n3 = 0;
            int n4 = indexIntersection.size();
            void v5 = var15_25;
            if(v5 == null)
            {
                throw new TypeCastException("null cannot be cast to non-null type java.lang.String");
            }
            Intrinsics.checkExpressionValueIsNotNull((Object) v5.substring(n3, n4), (String) "(this as java.lang.Strin\u2026ing(startIndex, endIndex)");
            collection.add(string);
        }
        Set rightCommon = CollectionsKt.toSet((Iterable) ((List) destination$iv$iv4));
        return Intrinsics.areEqual((Object) leftCommon, (Object) rightCommon);

         */
    }

    public static PhasedEvidence combineOverlapping(@NotNull PhasedEvidence left, @NotNull PhasedEvidence right)
    {
        return null;

        /*
        void $receiver$iv$iv;
        Intrinsics.checkParameterIsNotNull((Object) left, (String) "left");
        Intrinsics.checkParameterIsNotNull((Object) right, (String) "right");
        List indexUnion =
                CollectionsKt.sorted((Iterable) CollectionsKt.toList((Iterable) CollectionsKt.union((Iterable) ArraysKt.toSet((int[]) left.getAminoAcidIndices()), (Iterable) ArraysKt
                        .toSet((int[]) right.getAminoAcidIndices()))));
        List indexIntersection =
                CollectionsKt.sorted((Iterable) CollectionsKt.toList((Iterable) CollectionsKt.intersect((Iterable) ArraysKt.toSet((int[]) left
                        .getAminoAcidIndices()), (Iterable) ArraysKt.toSet((int[]) right.getAminoAcidIndices()))));
        int[] $receiver$iv = left.getAminoAcidIndices();
        Object object = $receiver$iv;
        Collection destination$iv$iv = new ArrayList();
        int n = ((void) $receiver$iv$iv).length;
        for(int i = 0; i < n; ++i)
        {
            void element$iv$iv;
            void it = element$iv$iv = $receiver$iv$iv[i];
            boolean bl = false;
            if(!(!indexIntersection.contains((int) it)))
            {
                continue;
            }
            destination$iv$iv.add((int) element$iv$iv);
        }
        List uniqueToLeft = (List) destination$iv$iv;
        Map result = new LinkedHashMap();
        Object object2 = left.getEvidence();
        Iterator<Map.Entry<String, Integer>> iterator = object2.entrySet().iterator();
        while(iterator.hasNext())
        {
            void leftSequence;
            Object element$iv$iv = object = (Object) iterator.next();
            object2 = (String) element$iv$iv.getKey();
            element$iv$iv = object;
            int leftCount = ((Number) element$iv$iv.getValue()).intValue();
            void it = leftSequence;
            int n2 = 0;
            int n3 = uniqueToLeft.size();
            void v0 = it;
            if(v0 == null)
            {
                throw new TypeCastException("null cannot be cast to non-null type java.lang.String");
            }
            Intrinsics.checkExpressionValueIsNotNull((Object) v0.substring(n2, n3), (String) "(this as java.lang.Strin\u2026ing(startIndex, endIndex)");
            Map.Entry<String, Integer> entry = leftSequence;
            n3 = uniqueToLeft.size();
            void v1 = entry;
            if(v1 == null)
            {
                throw new TypeCastException("null cannot be cast to non-null type java.lang.String");
            }
            Intrinsics.checkExpressionValueIsNotNull((Object) v1.substring(n3), (String) "(this as java.lang.String).substring(startIndex)");
            Object $i$f$filter = right.getEvidence();
            Iterator<Map.Entry<String, Integer>> iterator2 = $i$f$filter.entrySet().iterator();
            while(iterator2.hasNext())
            {
                String leftUnique;
                String rightCommon;
                String leftCommon;
                void rightSequence;
                Map.Entry<String, Integer> entry2 = entry = iterator2.next();
                $i$f$filter = entry2.getKey();
                entry2 = entry;
                int rightCount = ((Number) entry2.getValue()).intValue();
                void var18_24 = rightSequence;
                int n4 = 0;
                int n5 = indexIntersection.size();
                void v2 = var18_24;
                if(v2 == null)
                {
                    throw new TypeCastException("null cannot be cast to non-null type java.lang.String");
                }
                Intrinsics.checkExpressionValueIsNotNull((Object) v2.substring(n4, n5), (String) "(this as java.lang.Strin\u2026ing(startIndex, endIndex)");
                if(!Intrinsics.areEqual((Object) leftCommon, (Object) rightCommon))
                {
                    continue;
                }
                String combined = leftUnique + (String) rightSequence;
                result.put(combined, Math.min(leftCount, rightCount));
            }
        }
        return new PhasedEvidence(CollectionsKt.toIntArray((Collection) indexUnion), result);

         */
    }
}
