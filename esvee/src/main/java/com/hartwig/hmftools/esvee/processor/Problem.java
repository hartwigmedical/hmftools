package com.hartwig.hmftools.esvee.processor;

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
                Context == null ? "" : " while processing " + Context, Error);
    }
}
