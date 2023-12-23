package com.hartwig.hmftools.svassembly.assembly;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.svassembly.models.SupportedAssembly;
import com.hartwig.hmftools.svassembly.models.TrimmableAssembly;
import com.hartwig.hmftools.svassembly.util.Timeout;

import org.jetbrains.annotations.Nullable;

public enum AssemblyFiltering
{
    ;

    public static <T extends SupportedAssembly & TrimmableAssembly<T>> List<T> trimAndDeduplicate(final SupportChecker supportChecker,
            final List<T> assemblies, final Timeout timeout)
    {
        final List<T> trimmed = new ArrayList<>();
        assemblies.forEach(assembly -> trimmed.add(trimAssembly(assembly)));

        final Set<String> assemblyText = new HashSet<>();
        for(int i = 0; i < trimmed.size(); i++)
        {
            @Nullable
            final T assembly = trimmed.get(i);
            if(assembly != null && (!assemblyText.add(assembly.Assembly) || assembly.getSupportFragments().size() < 2))
                trimmed.set(i, null);
        }

        for (int i = 0; i < trimmed.size(); i++)
            for (int j = 0; j < trimmed.size(); j++)
            {
                if (j % 4 == 0)
                    timeout.checkTimeout();
                if (i == j)
                    continue;

                @Nullable
                final T left = trimmed.get(i);
                @Nullable
                final T right = trimmed.get(j);
                if (left == null || right == null)
                    continue;

                @Nullable
                final T deduped = resolveNearDuplicate(supportChecker, left, right);
                if (deduped == null)
                    continue;

                trimmed.set(i, deduped);
                trimmed.set(j, null);
            }

        return trimmed.stream()
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
    }

    @Nullable
    public static <T extends SupportedAssembly & TrimmableAssembly<T>> T trimAssembly(@Nullable final T assembly)
    {
        if(assembly == null)
            return null;
        else if(assembly.supportCount() == 0)
            return assembly;

        final var pair = assembly.computeBaseSupportAndContradiction();
        final int[] baseSupport = pair.getLeft();

        int rightTrimLength = 0;
        for(int i = assembly.Assembly.length() - 1; i >= 0; i--)
            if(baseSupport[i] == 0)
                rightTrimLength++;
            else
                break;

        int leftTrimLength = 0;
        for(int i = 0; i < assembly.Assembly.length(); i++)
            if(baseSupport[i] == 0)
                leftTrimLength++;
            else
                break;

        if(rightTrimLength <= 0 && leftTrimLength <= 0)
            return assembly;

        return assembly.trim(leftTrimLength, rightTrimLength);
    }

    @Nullable
    private static <T extends SupportedAssembly> T resolveNearDuplicate(final SupportChecker supportChecker, final T left, final T right)
    {
        // Work out if left/right are compatible
        @Nullable
        final Integer rightSupportingLeft = supportChecker.WeakSupport.supportIndex(left, right, 100);
        @Nullable
        final Integer leftSupportingRight = supportChecker.WeakSupport.supportIndex(right, left, 100);

        if(rightSupportingLeft == null && leftSupportingRight == null)
            return null;
        else if(rightSupportingLeft != null && leftSupportingRight == null)
            return left;
        else if(rightSupportingLeft == null) // && leftSupportingRight != null
            return right;
        assert rightSupportingLeft != null && leftSupportingRight != null;

        // TODO: Extend one assembly using the other

        if(left.getSupportFragments().size() > right.getSupportFragments().size())
            return left;
        else if(left.getSupportFragments().size() < right.getSupportFragments().size())
            return right;

        if(left.Assembly.length() > right.Assembly.length())
            return left;
        else if(left.Assembly.length() < right.Assembly.length())
            return right;

        return left; // Default to left, there is no reason for this.
    }
}
