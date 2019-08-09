package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicUser {

    @NotNull
    public abstract String username();

    @Nullable
    public abstract String areaOfExpertise();

    @NotNull
    public abstract CivicOrganization organization();

    @Nullable
    public abstract String twitterHandle();

    @NotNull
    public abstract String name();

    @Nullable
    public abstract String bio();

    @Nullable
    public abstract String url();

    @NotNull
    public abstract String createdAt();

    @NotNull
    public abstract CivicAvatars avatars();

    @Nullable
    public abstract String acceptedLicense();

    @Nullable
    public abstract String affiliation();

    @NotNull
    public abstract String avatarUrl();

    @NotNull
    public abstract String role();

    @Nullable
    public abstract String facebookProfile();

    @Nullable
    public abstract String linkedinProfile();

    @Nullable
    public abstract String orcid();

    @NotNull
    public abstract String displayName();

    @Nullable
    public abstract String lastSeenAt();

    @NotNull
    public abstract String featuredExpert();

    @NotNull
    public abstract String id();

    @Nullable
    public abstract String signupComplete();
}
