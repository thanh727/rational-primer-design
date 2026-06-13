const ms = parseInt(process.argv[2]) || 3000;
setTimeout(() => {
  process.exit(0);
}, ms);
