// 生成目录
function generateTOC() {
  const content = document.querySelector('.main-content');
  const headings = content.querySelectorAll('h1, h2, h3');
  const toc = document.querySelector('.table-of-contents ul');
  
  headings.forEach((heading, index) => {
    heading.id = `heading-${index}`;
    
    const li = document.createElement('li');
    const a = document.createElement('a');
    a.href = `#heading-${index}`;
    a.textContent = heading.textContent;
    a.style.paddingLeft = `${(heading.tagName[1] - 1) * 1}rem`;
    
    li.appendChild(a);
    toc.appendChild(li);
  });
}

// 控制目录的展开/收起
function initTOCBehavior() {
  const toc = document.querySelector('.table-of-contents');
  let isMouseInToc = false;

  toc.addEventListener('mouseenter', () => {
    isMouseInToc = true;
    toc.classList.add('expanded');
  });

  toc.addEventListener('mouseleave', (e) => {
    const tocRect = toc.getBoundingClientRect();
    const isInExtraArea = e.clientX >= tocRect.right && 
                         e.clientX <= tocRect.right + 20 && 
                         e.clientY >= tocRect.top && 
                         e.clientY <= tocRect.bottom;
    
    if (!isInExtraArea) {
      isMouseInToc = false;
      toc.classList.remove('expanded');
    }
  });
}

// 初始化
document.addEventListener('DOMContentLoaded', () => {
  generateTOC();
  initTOCBehavior();
}); 