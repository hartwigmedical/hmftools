package com.hartwig.hmftools.esvee.processor;

import com.hartwig.hmftools.esvee.JunctionProcessingException;

public class Problem
{
    public final String Message;
    public final Throwable Error;
    public final Object Context;

    public Problem(final String message, final Throwable error, final Object context)
    {
        Message = message;
        Error = error;
        Context = context;
    }

    @Override
    public String toString()
    {
        return String.format("%s%s: %s", Message,
                Context == null ? "" : " while processing " + Context,
                Error instanceof JunctionProcessingException ? Error.getClass() : Error);
    }
}
