
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - Gripss"

dependencies {
    implementation(project(":hmf-common"))
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.gripss.GripssApplication")
}
