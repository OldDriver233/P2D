<script setup lang="ts">
  import { ref } from 'vue'

  const count = ref(0)

</script>

<template>
  <div class="about">
    <h1>This is an about page</h1>
    <button @click="count++">{{ count }}</button>
  </div>
</template>

<style>
@media (min-width: 1024px) {
  .about {
    min-height: 100vh;
    display: flex;
    align-items: center;
  }
}
</style>
