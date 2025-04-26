// Ketcher 초기화
const ketcher = new window.Ketcher(document.getElementById('ketcher-container'));

// 이름 생성 함수
async function generateName() {
  const smiles = await ketcher.getSmiles();

  fetch('https://cactus.nci.nih.gov/chemical/structure/' + encodeURIComponent(smiles) + '/iupac_name')
    .then(response => response.text())
    .then(name => {
      document.getElementById('result').innerText = 'IUPAC 이름: ' + name;
    })
    .catch(error => {
      document.getElementById('result').innerText = '이름을 가져오는 데 실패했습니다.';
      console.error(error);
    });
}
