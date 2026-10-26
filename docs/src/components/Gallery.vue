<!-- The examples gallery: a grid of cards -->
<script setup lang="ts">
import GalleryImage, { type Props } from './GalleryImage.vue';

defineProps<{
    images: Props[];
}>();
</script>

<template>
    <div class="gallery">
        <GalleryImage v-for="image in images" v-bind="image" />
    </div>
</template>

<style scoped>
.gallery {
    display: flex;
    flex-wrap: wrap;
    justify-content: flex-start;
    margin: 0 -12px;
}
</style>
